#include "Glial.h"
#include "Eigen/Dense"
#include <Eigen/Geometry>
#include <Eigen/Core>
#include "constants.h"
#include <numeric>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <set>
#include <limits>
#include <cmath>
#include <limits>
#include <algorithm> // for std::min initializer_list


using namespace Eigen;
using namespace std;

Glial::Glial()
{}

Glial::~Glial()
{}

Glial::Glial(const Glial &gl)
{
    id = gl.id;
    soma = gl.soma;
    processes = gl.processes;
    grid = gl.grid;
    percolation = gl.percolation;
    prob_cross_e_i = gl.prob_cross_e_i;
    prob_cross_i_e = gl.prob_cross_i_e;
    diffusivity_i = gl.diffusivity_i;
    diffusivity_e = gl.diffusivity_e;
    count_perc_crossings = gl.count_perc_crossings;

};


inline bool Glial::is_empty(const Box& b) {
    return b.x_min > b.x_max || b.y_min > b.y_max || b.z_min > b.z_max;
}
inline void Glial::extend(Box& b, const Eigen::Vector3d& p) {
    if (is_empty(b)) { b = {p.x(),p.x(),p.y(),p.y(),p.z(),p.z()}; return; }
    b.x_min = std::min(b.x_min, p.x()); b.x_max = std::max(b.x_max, p.x());
    b.y_min = std::min(b.y_min, p.y()); b.y_max = std::max(b.y_max, p.y());
    b.z_min = std::min(b.z_min, p.z()); b.z_max = std::max(b.z_max, p.z());
}
// 3-int hash for buckets
inline uint64_t hash3(int x, int y, int z) {
    // simple 64-bit mix (SplitMix-like)
    auto mix = [](uint64_t v){ v += 0x9e3779b97f4a7c15ULL; v = (v^(v>>30))*0xbf58476d1ce4e5b9ULL;
                               v = (v^(v>>27))*0x94d049bb133111ebULL; return v^(v>>31); };
    return mix((uint64_t)(uint32_t)x) ^ (mix((uint64_t)(uint32_t)y)<<1) ^ (mix((uint64_t)(uint32_t)z)<<2);
}
void Glial::build_glia_grid_processes(const std::vector<std::vector<Sphere>>& processes,
                                      double cell_size, double pad)
{
    // init grid
    grid = HashGrid{};
    grid.cell      = (cell_size > 0.0 ? cell_size : 1.0);
    grid.build_pad = std::max(0.0, pad);              // remember build pad
    grid.max_radius_plus_pad = 0.0;                   // start tracking

    // 1) compute big box from processes, already padded by `pad`
    Box B = make_empty_box();
    for (const auto& branch : processes) {
        for (const auto& s : branch) {
            if (s.radius <= 0.0) continue;
            const double R = s.radius + pad;
            extend(B, s.P - Eigen::Vector3d::Constant(R));
            extend(B, s.P + Eigen::Vector3d::Constant(R));
            // track global max radius+pad
            if (R > grid.max_radius_plus_pad) grid.max_radius_plus_pad = R;
        }
    }

    // Also extend by soma (helpful if processes are sparse/empty)
    {
        const double R = soma.radius + pad;
        extend(B, soma.P - Eigen::Vector3d::Constant(R));
        extend(B, soma.P + Eigen::Vector3d::Constant(R));
        if (R > grid.max_radius_plus_pad) grid.max_radius_plus_pad = R;
    }

    // Fallback if still empty (degenerate case)
    if (is_empty(B)) {
        const double R = std::max(1e-6, soma.radius + pad);
        const Eigen::Vector3d mn = soma.P - Eigen::Vector3d::Constant(R);
        const Eigen::Vector3d mx = soma.P + Eigen::Vector3d::Constant(R);
        B = {mn.x(), mx.x(), mn.y(), mx.y(), mn.z(), mx.z()};
        grid.max_radius_plus_pad = std::max(grid.max_radius_plus_pad, R);
    }

    // 2) set grid origin and store big box
    grid.origin  = Eigen::Vector3d(B.x_min, B.y_min, B.z_min);
    grid.big_box = B;

    // (Optional) reserve to avoid reallocation
    size_t total = 0;
    for (const auto& br : processes) total += br.size();
    grid.objs.reserve(total);

    // 3) add spheres into buckets — store indices (b,i)
    auto add = [&](int b, int i) {
        const Sphere& s = processes[b][i];
        if (s.radius <= 0.0) return;
        const double R = s.radius + pad;
        const Eigen::Vector3d mn = s.P - Eigen::Vector3d::Constant(R);
        const Eigen::Vector3d mx = s.P + Eigen::Vector3d::Constant(R);

        const Eigen::Array3i imin = ((mn - grid.origin).array() / grid.cell).floor().cast<int>();
        const Eigen::Array3i imax = ((mx - grid.origin).array() / grid.cell).floor().cast<int>();

        const int idx = static_cast<int>(grid.objs.size());
        grid.objs.push_back({b, i});     // store reference into `processes`

        for (int ix = imin.x(); ix <= imax.x(); ++ix)
          for (int iy = imin.y(); iy <= imax.y(); ++iy)
            for (int iz = imin.z(); iz <= imax.z(); ++iz)
              grid.buckets[hash3(ix,iy,iz)].push_back(idx);
    };

    for (int b = 0; b < (int)processes.size(); ++b) {
        const auto& branch = processes[b];
        for (int i = 0; i < (int)branch.size(); ++i) {
            const Sphere& s = branch[i];
            if (s.radius <= 0.0) continue;
            add(b, i);
        }
    }

    // 4) finalize derived value
    grid.max_r_cells = std::max(
        1,
        (int)std::ceil(grid.max_radius_plus_pad / grid.cell)
    );
}



inline bool Glial::point_in_inflated_aabb(const Eigen::Vector3d& p,
                                   double d)
{
    Box box = grid.big_box;
    const double infl = std::max(0.0, d);
    return (p.x() >= box.x_min - infl && p.x() <= box.x_max + infl &&
            p.y() >= box.y_min - infl && p.y() <= box.y_max + infl &&
            p.z() >= box.z_min - infl && p.z() <= box.z_max + infl);
}

void Glial::set_up_glialcell(std::vector<Sphere> &spheres_to_add) {

    // Clear existing boxes and initialize variables
    processes.clear();

    if (spheres_to_add.empty()) {
        return;
    }

    for (const auto &sphere : spheres_to_add) {

        int branch_id = sphere.branch_id;
        int cell_id = sphere.id;

        if (cell_id < 0|| branch_id < 0) {
            cerr << "Error: Sphere with invalid id or branch_id: " << cell_id << ", " << branch_id << endl;
            assert(0);
        }

        if (branch_id >= processes.size()) {
            processes.resize(branch_id + 1);
        }
        processes[branch_id].push_back(sphere);
        
    }
    for (auto &branch : processes) {
        for (auto &sphere : branch) {
            if (sphere.id < 0 || sphere.branch_id < 0) {
                cerr << "Error in build: Sphere with invalid id or branch_id: " << sphere.id << ", " << sphere.branch_id << endl;
                assert(0);
            }
        }
    }
    double grid_cell_size = 0.0005;
    double pad = barrier_tickness;
    build_glia_grid_processes(processes,grid_cell_size, pad);

}

bool Glial::is_point_near_glia(const Eigen::Vector3d& p,
                        double d)
{
    if ((p - soma.P).squaredNorm() <= (soma.radius + d)*(soma.radius + d))
        return true;
    if (!point_in_inflated_aabb(p, d)) return false;

    Eigen::Array3i ic = ((p - grid.origin).array() / grid.cell).floor().cast<int>();
    const int L = std::max(1, (int)std::ceil(d / grid.cell));
    const double d2 = d*d;

    for (int dx=-L; dx<=L; ++dx)
      for (int dy=-L; dy<=L; ++dy)
        for (int dz=-L; dz<=L; ++dz) {
          auto it = grid.buckets.find(hash3(ic[0]+dx, ic[1]+dy, ic[2]+dz));
          if (it == grid.buckets.end()) continue;
          for (int idx : it->second) {
            auto [b, i] = grid.objs[idx];  // (branch_id, cell_id)
            const Sphere* s = &processes[b][i];
            const Eigen::Vector3d v = p - s->P;
            if (v.squaredNorm() <= (s->radius + d)*(s->radius + d)) return true;
          }
        }
    return false;
}


// Optional: segment vs AABB to clamp traversal to the grid big_box (slab method)
inline bool Glial::segment_aabb_intersect(const Eigen::Vector3d& p0,
                                   const Eigen::Vector3d& p1,
                                   const Box& box,
                                   double& tEnter, double& tExit)
{
    const Eigen::Vector3d d = p1 - p0;      // segment direction
    tEnter = 0.0;                            // param in [0,1]
    tExit  = 1.0;

    for (int a = 0; a < 3; ++a) {
        double minA;
        double maxA;

        if (a==0){
            minA = box.x_min; maxA = box.x_max;
        } else if (a==1) {
            minA = box.y_min; maxA = box.y_max;
        } else { // a==2
            minA = box.z_min; maxA = box.z_max;
        }

        if (std::abs(d[a]) < 1e-15) {
            // Segment is parallel to this axis' slabs: must be within slab
            if (p0[a] < minA || p0[a] > maxA) return false;
            continue;
        }

        const double inv = 1.0 / d[a];
        double t0 = (minA - p0[a]) * inv;
        double t1 = (maxA - p0[a]) * inv;
        if (t0 > t1) std::swap(t0, t1);

        tEnter = std::max(tEnter, t0);
        tExit  = std::min(tExit,  t1);
        if (tEnter > tExit) return false;   // no overlap
    }
    return true;  // segment intersects the box for t in [tEnter, tExit]
}

// Gathers process-sphere candidate indices from the grid cells crossed by the step.
// - p0: segment start
// - dir_unit: unit direction
// - L: segment length (|step|)
// - out_ids: candidate indices (into G.objs), deduplicated
void Glial::gather_candidates_DDA(const Eigen::Vector3d& p0,
                           const Eigen::Vector3d& dir_unit,
                           double L,
                           std::vector<int>& out_ids)
{
    out_ids.clear();
    if (L <= 0) return;

    const double eps = 1e-12;

    // Clamp traversal to grid big_box (avoid walking forever if outside)
    double tEnter=0.0, tExit=1.0;
    const Eigen::Vector3d p1 = p0 + dir_unit * L;
    if (!segment_aabb_intersect(p0, p1, grid.big_box, tEnter, tExit)) {
        return; // segment misses overall grid AABB
    }

    // Convert tEnter/tExit in [0,1] to [0,L]
    double t0 = std::max(0.0, tEnter * L);
    double t1 = std::min(L,     tExit  * L);
    if (t0 > t1) return;

    // Start position (nudged slightly inside the grid cell to avoid exact-boundary issues)
    Eigen::Vector3d p = p0 + dir_unit * (t0 + eps);

    // Current integer cell
    auto to_cell = [&](const Eigen::Vector3d& q)->Eigen::Array3i {
        return ((q - grid.origin).array() / grid.cell).floor().cast<int>();
    };
    Eigen::Array3i cell = to_cell(p); // current cell

    // Step direction per axis
    Eigen::Array3i step;
    for (int a=0; a<3; ++a)
        step[a] = (dir_unit[a] > 0) ? 1 : (dir_unit[a] < 0 ? -1 : 0);

    // Min corner of this cell
    Eigen::Array3d cellMin = (grid.origin.array() + cell.cast<double>() * grid.cell);

    // Compute tMax (param distance to next boundary) for each axis
    Eigen::Array3d tMax;
    for (int a=0; a<3; ++a) {
        if (step[a] > 0) {
            double nextFace = cellMin[a] + grid.cell;
            tMax[a] = (nextFace - p[a]) / (dir_unit[a] + ((dir_unit[a]==0)?eps:0));
        } else if (step[a] < 0) {
            double prevFace = cellMin[a];
            tMax[a] = (prevFace - p[a]) / (dir_unit[a] + ((dir_unit[a]==0)?eps:0));
        } else {
            tMax[a] = std::numeric_limits<double>::infinity();
        }
        if (tMax[a] < 0) tMax[a] = 0; // guard tiny negatives due to eps nudge
    }

    // Distance in t to cross one full cell on each axis
    Eigen::Array3d tDelta;
    for (int a=0; a<3; ++a) {
        if (step[a] != 0)
            tDelta[a] = grid.cell / std::abs(dir_unit[a]);
        else
            tDelta[a] = std::numeric_limits<double>::infinity();
    }

    // Traverse
    std::unordered_set<int> seen; seen.reserve(128);
    double t = t0;

    while (t <= t1 + eps) {
        // 1) Gather sphere indices from current cell
        auto it = grid.buckets.find(hash3(cell[0], cell[1], cell[2]));
        if (it != grid.buckets.end()) {
            for (int idx : it->second)
                if (seen.insert(idx).second) out_ids.push_back(idx);
        }

        // 2) Advance to next cell boundary (smallest tMax)
        int axis = 0;
        if (tMax[1] < tMax[axis]) axis = 1;
        if (tMax[2] < tMax[axis]) axis = 2;

        double tNext = t + tMax[axis];
        if (tNext > t1 + eps) break;  // next boundary beyond segment end

        // Step to next cell on that axis
        cell[axis] += step[axis];
        // Shift reference point and tMax: after crossing, the next boundary along this axis is tMax += tDelta
        t     = tNext;
        tMax -= Eigen::Array3d::Constant(tMax[axis]); // zero the chosen axis tMax
        tMax[axis] += tDelta[axis];

        // (Optional) early-out if cell is far outside big_box; with slab clamp above this is rare.
    }
}



// Return true if the segment p0 + t*dir_unit, t∈(0,L] intersects a sphere (C,R).
// t_enter <= t_exit are clamped to [0, L].
inline bool Glial::raySphere(const Eigen::Vector3d& p0,
                      const Eigen::Vector3d& dir_unit, // must be unit
                      const Eigen::Vector3d& C,
                      double R,
                      double& t_enter,
                      double& t_exit)
{
    // Solve ||(p0 - C) + t*dir||^2 = R^2  with a=1 (dir is unit)
    const Eigen::Vector3d oc = p0 - C;
    const double b = oc.dot(dir_unit);
    const double c = oc.squaredNorm() - R*R;
    const double disc = b*b - c;
    if (disc < 0.0) {
        //cout << "No intersection with sphere: disc=" << disc << endl;
        return false;
    }

    const double s = std::sqrt(std::max(0.0, disc));
    double t0 = -b - s;    // enter
    double t1 = -b + s;    // exit
    if (t0 > t1) std::swap(t0, t1);
    t_enter = t0;
    t_exit  = t1;
    return true;
}

bool Glial::checkCollision(const Walker& walker,
                           Eigen::Vector3d& step,
                           const double& step_length,
                           Collision& collision)
{

    // Normalize direction
    const double L = (step_length > 0.0) ? step_length : step.norm();
    if (L <= 0.0) { 
        collision.type = Collision::null; 
        return false; 
        }
    const Eigen::Vector3d dir = step.normalized();
    const Eigen::Vector3d p0  = walker.pos_v;


    const double Rpad = grid.build_pad;
    int occ0 = occupancy_at_point(p0, Rpad);
    bool is_inside = (occ0 > 0);
    const bool start_inside = (walker.location == Walker::intra);

    if (start_inside && occ0 == 0) {
        const double bias = std::max(1e-6, 1e-3 * grid.cell);
        double sd = signed_distance_to_union(p0, Rpad); // negative means inside
        if (sd < bias) {
            // Treat as inside for this step; you are glued to the surface
            occ0 = 1;
            std::cerr << "Sticky-inside: sd=" << sd << "\n";
        }
    }
    if (start_inside && !is_inside){
        // problem
        collision.type = Collision::hit;
        collision.col_location  = Collision::outside;
        collision.perm_crossing = 0.0;
        cout << "Error: walker started inside glia but is not inside any sphere at p0\n";
        cout << "p0=" << p0.transpose() << " dir=" << dir.transpose() 
             << " occ0=" << occ0 <<  "\n";
        //assert(0);
        return true;
    }
    else if (!start_inside && is_inside){
        // problem
        collision.type = Collision::hit;
        collision.col_location  = Collision::inside;
        collision.perm_crossing = 0.0;
        cout << "Error: walker started outside glia but is inside a sphere at p0\n";
        return true;
    }

    //cout <<"----- Walker checkCollision at p0=" << p0.transpose() << " dir=" << dir.transpose() << " L=" << L << "\n";
    /*
    bool inside_glial = isPosInsideGlialCell(p0, 1e-3);
    if (!inside_glial) {
        cout <<"Walker is outside glia at p0=" << p0.transpose() << "\n";
        assert(0);
    }
    */
    

    // Gather candidates
    std::vector<int> cand_ids;
    gather_candidates_DDA(p0, dir, L, cand_ids);

    struct Ev { double t; int delta; const Sphere* s; };
    std::vector<Ev> evs; evs.reserve(cand_ids.size()*2 + 2);

    const double epsT = std::max(1e-12, 1e-6 * grid.cell);

    auto addSphereEvents = [&](const Sphere* s,
                            const Eigen::Vector3d& p0,
                            const Eigen::Vector3d& dir) {
        if (!s) return;

        double t0, t1;
        const double Rin = s->radius + Rpad;
        if (!raySphere(p0, dir, s->P, Rin, t0, t1)) return;

        // Ensure t0 <= t1 (if your raySphere doesn’t guarantee it)
        if (t1 < t0) std::swap(t0, t1);

        // Completely outside the segment?
        if (t1 < 0.0 || t0 > L) return;

        // Clip to [0, L]
        t0 = std::max(0.0, std::min(t0, L));
        t1 = std::max(0.0, std::min(t1, L));
        if (t1 <= t0 + epsT) return;

        const bool inside0 = (p0 - s->P).squaredNorm() <= Rin*Rin + 1e-12;

        if (!inside0) {
            if (t0 > epsT) evs.push_back({t0, +1, s});  // ENTER
            if (t1 > epsT) evs.push_back({t1, -1, s});  // EXIT
        } else {
            if (t1 > epsT) evs.push_back({t1, -1, s});  // EXIT only
        }
    };


    addSphereEvents(&soma, p0, dir);
    //cout <<"cand_ids.size()=" << cand_ids.size() << " grid.objs.size()=" << grid.objs.size() << "\n";

    // SAFETY: validate every candidate index + pointer
    size_t bad_idx = 0, null_ptr = 0;
    for (size_t k = 0; k < cand_ids.size(); ++k) {
        int idx = cand_ids[k];
        if (idx < 0 || static_cast<size_t>(idx) >= grid.objs.size()) {
            ++bad_idx;
            std::cerr << "BAD candidate idx=" << idx
                      << " (objs.size()=" << grid.objs.size()
                      << ", k=" << k << ")\n";
            continue;
        }
        auto [b, i] = grid.objs[idx];   
        const Sphere* sp = &processes[b][i];
        if (!sp) { ++null_ptr; std::cerr << "NULL grid.objs["<<idx<<"]\n"; continue; }
        addSphereEvents(sp, p0, dir);
    }
    if (bad_idx || null_ptr) {
        std::cerr << "Summary: bad_idx=" << bad_idx << " null_ptr=" << null_ptr << "\n";
    }

    if (evs.empty()) { 
        collision.type = Collision::null; 
        //cout << "No sphere collisions found for walker at p0=" << p0.transpose() << " with direction :"<< dir.transpose() <<" and L : "<< L << "\n";
        return false; 
    }

    std::sort(evs.begin(), evs.end(), [](const Ev& a, const Ev& b){ return a.t < b.t; });

    // Drop any events beyond L (in case raySphere or FP noise sneaks one in)
    while (!evs.empty() && evs.back().t > L + 1e-12) evs.pop_back();
    if (evs.empty()) { collision.type = Collision::null; return false; }

    // 5) Sweep to find first union boundary:
    
    int occ = occ0;
    

    const Ev* hit = nullptr;

    if (!start_inside) {
        for (const auto& e : evs) { occ += e.delta; if (occ > 0) { hit = &e; break; } }
    } else {
        for (const auto& e : evs) { occ += e.delta; if (occ == 0) { hit = &e; break; } }
    }

    if (!hit) { 
        collision.type = Collision::null; 
        //cout << "No sphere collisions found for walker at p0=" << p0.transpose() 
        //     << " with occupancy start=" << occ0 << " final=" << occ << "\n";
        return false; 
    }

    // 6) Fill collision info
    double t_hit = hit->t - Rpad;

    const Eigen::Vector3d pos = p0 + t_hit * dir;
    // is the hit point inside the sphere?
    bool next_pos_inside = isPosInsideGlialCell(pos, 0);
    if (!next_pos_inside && start_inside){
        double old_t_hit = t_hit;
        bool ok = ensure_same_compartment_at_hit(p0, dir, start_inside, Rpad, grid.cell, t_hit);
        if (!ok) {
            // We were already on the wrong side at t=0; do a t=0 bounce
            // (or your preferred “repair” behavior)
            t_hit = 0.0;
            cout << "Error: ensure_same_compartment_at_hit failed for walker at p0=" 
                 << p0.transpose() << " dir=" << dir.transpose() 
                 << " start_inside=" << start_inside << " occ0=" << occ0 
                 << " occ_final=" << occ << "\n";
            assert(0);
        }
        /*
        cout << "Warning: ensure_same_compartment_at_hit repaired t_hit=" 
             << t_hit << ", old t_hit : "<< old_t_hit <<" for walker at p0=" << p0.transpose() 
             << " dir=" << dir.transpose() << " start_inside=" << start_inside 
             << " occ0=" << occ0 << " occ_final=" << occ << "\n";
             */
    }
    const Eigen::Vector3d n   = (pos - hit->s->P).normalized();
    const double dn = dir.dot(n);
    const Eigen::Vector3d bounced = dir - 2.0 * dn * n;

    collision.type = Collision::hit;
    collision.collision_point = pos;
    collision.t = t_hit;
    collision.bounced_direction = bounced.normalized();
    collision.obstacle_type = 1;                    // your code’s convention
    collision.obstacle_ind  = hit->s->id;           // or branch/local as needed
    collision.col_location  = start_inside ? Collision::inside : Collision::outside;
    collision.perm_crossing = 0.0;
    /*
    cout <<"hit sphere id=" << collision.obstacle_ind
         << " at t=" << t_hit << " pos=" << pos.transpose()
         << " dir=" << dir.transpose() << " bounced=" << bounced.transpose()
         << " col_loc=" << collision.col_location
         << " start_inside=" << start_inside
         << " occ0=" << occ0
         << " occ_final=" << occ
         << "\n";
         */
         

    // Optional permeability
    if (percolation > 0.0) {
        static thread_local std::mt19937 gen{std::random_device{}()};
        std::uniform_real_distribution<double> U(0.0,1.0);
        const double p_cross = start_inside ? prob_cross_i_e : prob_cross_e_i;
        if (U(gen) < p_cross) {
            collision.perm_crossing = p_cross;
            collision.bounced_direction = dir; // continue forward
        }
    }
    return true;
}

// Returns true if we found a t' <= t_hit such that occupancy(p0+dir*t') == start_inside.
// t_hit is modified in place. Uses bisection with a tolerance tied to the grid cell.
bool Glial::ensure_same_compartment_at_hit(const Eigen::Vector3d& p0,
                                           const Eigen::Vector3d& dir_unit,
                                           bool start_inside,
                                           double pad,      // use grid.build_pad
                                           double cell,     // grid.cell
                                           double& t_hit,   // in/out
                                           int max_iter) 
{
    const double tol = std::max(1e-12, 1e-6 * cell);

    auto occ = [&](double t)->bool {
        return occupancy_at_point(p0 + dir_unit * t, pad) > 0;
    };

    // If just stepping back a hair already fixes it, do that quickly.
    double t_try = std::max(0.0, t_hit - tol);
    if (occ(t_try) == start_inside) { t_hit = t_try; return true; }

    // Otherwise, bisection on [0, t_hit] to find largest t with desired occupancy.
    if (occ(0.0) != start_inside) {
        // This means we started on the wrong side; caller can handle as "t=0 bounce".
        t_hit = 0.0;
        return false;
    }

    double lo = 0.0, hi = t_hit, best = 0.0;
    for (int it = 0; it < max_iter; ++it) {
        double mid = 0.5 * (lo + hi);
        if (occ(mid) == start_inside) { best = mid; lo = mid; }
        else                          { hi   = mid; }
        if (hi - lo <= tol) break;
    }
    t_hit = std::max(0.0, best - 0.5 * tol); // nudge a touch inward for safety
    return true;
}

int Glial::occupancy_at_point(const Eigen::Vector3d& p, double margin) const
{
    const double pad = std::max(0.0, margin);
    int occ = 0;

    // Soma
    {
        const double R2 = (soma.radius + pad) * (soma.radius + pad);
        if ((p - soma.P).squaredNorm() <= R2 + 1e-12) ++occ;
    }

    // Big-box quick reject
    const Box& B = grid.big_box;
    if (B.x_min > B.x_max || B.y_min > B.y_max || B.z_min > B.z_max) return occ;
    if (p.x() < B.x_min - pad || p.x() > B.x_max + pad ||
        p.y() < B.y_min - pad || p.y() > B.y_max + pad ||
        p.z() < B.z_min - pad || p.z() > B.z_max + pad) return occ;

    const double cell = grid.cell;
    const Eigen::Array3i ic = ((p - grid.origin).array() / cell).floor().cast<int>();

    auto scanNeighbors = [&](int Lc) {
        std::unordered_set<int> seen; seen.reserve(64);
        for (int dx=-Lc; dx<=Lc; ++dx)
          for (int dy=-Lc; dy<=Lc; ++dy)
            for (int dz=-Lc; dz<=Lc; ++dz) {
                auto it = grid.buckets.find(hash3(ic[0]+dx, ic[1]+dy, ic[2]+dz)); // find cell in neighborhood
                if (it == grid.buckets.end()) continue; // if cell isnt within any sphere continue
                for (int idx : it->second) { // for each sphere in the cell
                    if (!seen.insert(idx).second) continue; // deduplicate
                    auto [b,i] = grid.objs[idx]; // find sphere
                    if ((size_t)b >= processes.size() || (size_t)i >= processes[b].size()) continue;
                    const Sphere& s = processes[b][i]; // find sphere
                    const double R2 = (s.radius + pad) * (s.radius + pad);
                    if ((p - s.P).squaredNorm() <= R2 + 1e-12) ++occ;
                }
            }
    };

    // Base radius in cells: cover pad + tiny safety
    int Lc = std::max(1, (int)std::ceil((pad + 1e-3*cell) / cell));
    scanNeighbors(Lc);

    // If nothing found, escalate once; this helps when p is right on a cell boundary
    if (occ == 0) scanNeighbors(Lc + 1);

    return occ;
}

double Glial::signed_distance_to_union(const Eigen::Vector3d& p, double margin) const
{
    const double pad = std::max(0.0, margin);
    double dmin = std::numeric_limits<double>::infinity();

    auto upd = [&](const Sphere& s){
        const double d = (p - s.P).norm() - (s.radius + pad);
        if (d < dmin) dmin = d;
    };

    // Soma
    upd(soma);

    // Big-box quick reject: if outside by more than dmin, you can early-return,
    // but we keep it simple here and just scan neighbors like occupancy:
    const double cell = grid.cell;
    const Eigen::Array3i ic = ((p - grid.origin).array() / cell).floor().cast<int>();
    const int Lc = std::max(1, (int)std::ceil((pad + 1e-3*cell) / cell));

    std::unordered_set<int> seen; seen.reserve(64);
    for (int dx=-Lc; dx<=Lc; ++dx)
      for (int dy=-Lc; dy<=Lc; ++dy)
        for (int dz=-Lc; dz<=Lc; ++dz) {
          auto it = grid.buckets.find(hash3(ic[0]+dx, ic[1]+dy, ic[2]+dz));
          if (it == grid.buckets.end()) continue;
          for (int idx : it->second) {
            if (!seen.insert(idx).second) continue;
            auto [b,i] = grid.objs[idx];
            if ((size_t)b >= processes.size() || (size_t)i >= processes[b].size()) continue;
            upd(processes[b][i]);
          }
        }
    return dmin; // <0 inside, >0 outside, ~0 on surface (w.r.t. margin)
}

void Glial::set_prob_crossings(double step_length_pref){

    double prob_cross_i_e_, prob_cross_e_i_;
    double dse, dsi;

    if (percolation > 0.0){
        
        dse = sqrt(step_length_pref*this->diffusivity_e);
        dsi = sqrt(step_length_pref*this->diffusivity_i);

        prob_cross_i_e_ = percolation * dsi * 2. / 3. / this->diffusivity_i;
        prob_cross_e_i_ = percolation * dse * 2. / 3. / this->diffusivity_e; 

        this->prob_cross_e_i = prob_cross_e_i_ / (1.+ 0.5 * (prob_cross_e_i_ + prob_cross_i_e_));
        this->prob_cross_i_e = prob_cross_i_e_ / (1.+ 0.5 * (prob_cross_e_i_ + prob_cross_i_e_));
    }
            
    // for all spheres
    for(unsigned i= 0 ; i < processes.size();i++){
        for (unsigned j = 0; j < processes[i].size(); j++){
            processes[i][j].percolation = this->percolation;
            if(processes[i][j].percolation > 0.0){

                processes[i][j].prob_cross_e_i = this->prob_cross_e_i;
                processes[i][j].prob_cross_i_e = this->prob_cross_i_e;

                processes[i][j].diffusivity_e = this->diffusivity_e;
                processes[i][j].diffusivity_i = this->diffusivity_i;
            }
        }
    }

    // soma 
    soma.percolation = this->percolation;
    if(soma.percolation > 0.0){
        soma.prob_cross_e_i = this->prob_cross_e_i;
        soma.prob_cross_i_e = this->prob_cross_i_e;

        soma.diffusivity_e = this->diffusivity_e;
        soma.diffusivity_i = this->diffusivity_i;
        
    }
}



double Glial::minDistance(const Walker& w) const
{
    const Box& box = grid.big_box;
    // if you have a helper is_empty(box), use that; otherwise:
    if (box.x_min > box.x_max || box.y_min > box.y_max || box.z_min > box.z_max)
        return std::numeric_limits<double>::infinity();

    const Eigen::Vector3d& p = w.pos_v;

    // Closest point on the box to p (clamp)
    Eigen::Vector3d q;
    q.x() = std::min(std::max(p.x(), box.x_min), box.x_max);
    q.y() = std::min(std::max(p.y(), box.y_min), box.y_max);
    q.z() = std::min(std::max(p.z(), box.z_min), box.z_max);

    const double d2 = (p - q).squaredNorm();
    if (d2 > 0.0) {
        // Outside: Euclidean distance to the box
        return std::sqrt(d2);
    } else {
        // Inside: minimum distance to any face (depth to surface)
        const double dx_min = p.x() - box.x_min;
        const double dx_max = box.x_max - p.x();
        const double dy_min = p.y() - box.y_min;
        const double dy_max = box.y_max - p.y();
        const double dz_min = p.z() - box.z_min;
        const double dz_max = box.z_max - p.z();
        return std::min({dx_min, dx_max, dy_min, dy_max, dz_min, dz_max});
    }
}

bool Glial::isPosInsideGlialCell(const Eigen::Vector3d& p, double margin) const
{
    const double pad = margin;

    // 0) Soma
    {
        const double R = soma.radius + pad;
        if ((p - soma.P).squaredNorm() <= R*R) return true;
    }

    // 1) Big-box quick reject
    const Box& B = grid.big_box;
    if (B.x_min > B.x_max || B.y_min > B.y_max || B.z_min > B.z_max) return false;
    if (p.x() < B.x_min - pad || p.x() > B.x_max + pad ||
        p.y() < B.y_min - pad || p.y() > B.y_max + pad ||
        p.z() < B.z_min - pad || p.z() > B.z_max + pad) return false;

    // 2) Cell index
    const Eigen::Array3i ic = ((p - grid.origin).array() / grid.cell).floor().cast<int>();

    auto scanNeighbors = [&](int Lc)->bool {
        std::unordered_set<int> seen; seen.reserve(64);
        for (int dx = -Lc; dx <= Lc; ++dx)
          for (int dy = -Lc; dy <= Lc; ++dy)
            for (int dz = -Lc; dz <= Lc; ++dz) {
                const uint64_t key = hash3(ic[0]+dx, ic[1]+dy, ic[2]+dz);
                auto it = grid.buckets.find(key);
                if (it == grid.buckets.end()) continue;
                for (int idx : it->second) {
                    if (!seen.insert(idx).second) continue;
                    auto [b, i] = grid.objs[idx];
                    if ((size_t)b >= processes.size() || (size_t)i >= processes[b].size()) continue;
                    const Sphere& s = processes[b][i];
                    const double R = s.radius + pad;
                    if ((p - s.P).squaredNorm() <= R*R) return true;
                }
            }
        return false;
    };

    // Neighborhood radius in cells: cover pad + a tiny safety
    int Lc = std::max(1, (int)std::ceil((pad + 1e-3*grid.cell) / grid.cell));
    if (scanNeighbors(Lc)) return true;

    // Fallback: escalate once (helps when grid.cell is *very* small)
    if (scanNeighbors(Lc+1)) return true;

    return false;
}
