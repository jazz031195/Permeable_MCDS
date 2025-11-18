
#include "Blood_Vessel.h"
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

Blood_Vessel::Blood_Vessel()
{}

Blood_Vessel::~Blood_Vessel()
{}

Blood_Vessel::Blood_Vessel(const Blood_Vessel &bv)
{
    id = bv.id;
    radius = bv.radius;
    pressure_diff = bv.pressure_diff;
    viscosity = bv.viscosity;
    flow = bv.flow;
    max_velocity = bv.max_velocity;
    spheres = bv.spheres;
    skeleton = bv.skeleton;
    grid = bv.grid;

};

void Blood_Vessel::distance_to_skeleton(const Walker &w, double& min_dist, Eigen::Vector3d& tangent){

    // distance O to the vessel skeleton
    Eigen::Vector3d O = w.pos_v;
    min_dist = std::numeric_limits<double>::max();
    int index = -1;
    for (unsigned i = 0; i < skeleton.size(); i++){
        double dist = (O - skeleton[i]).norm();
        if (dist < min_dist){
            min_dist = dist;
            index = i;
        }
    }


    if (w.normal[2] == 0){ // before colliding with voxel wall

        // tangent direction at the closest point
        if (index == 0){
            tangent = (skeleton[1] - skeleton[0]).normalized();
        }
        else if (index == skeleton.size() - 1){
            tangent = (skeleton[skeleton.size() - 1] - skeleton[skeleton.size() - 2]).normalized();
        }
        else{
            tangent = (skeleton[index + 1] - skeleton[index - 1]).normalized();
        }
    }
    else{ // after colliding with voxel wall
        // tangent direction at the closest point
        if (index == 0){
            tangent = (skeleton[0] - skeleton[1]).normalized();
        }
        else if (index == skeleton.size() - 1){
            tangent = (skeleton[skeleton.size() - 2] - skeleton[skeleton.size() - 1]).normalized();
        }
        else{
            tangent = (skeleton[index - 1] - skeleton[index + 1]).normalized();
        }
    }
}

void Blood_Vessel::velocity(const Walker &w, double& v, Eigen::Vector3d& flow_direction){

    double min_dist;

    Eigen::Vector3d tangent;

    distance_to_skeleton(w, min_dist, tangent);


    if (min_dist >= this->radius) {
        v = 0.0;
        flow_direction = Eigen::Vector3d::Zero();
        return;
    }

    double radial_distance = min_dist; // distance is negative inside the vessel

    v= (pressure_diff/(4*viscosity))*(this->radius*this->radius - radial_distance*radial_distance); // Parabolic profile

    flow_direction = tangent;

}



inline bool Blood_Vessel::is_empty(const Box& b) {
    return b.x_min > b.x_max || b.y_min > b.y_max || b.z_min > b.z_max;
}
inline void Blood_Vessel::extend(Box& b, const Eigen::Vector3d& p) {
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


void Blood_Vessel::set_spheres(std::vector<Sphere> &spheres_to_add) {

    // Clear existing boxes and initialize variables
    spheres.clear();

    if (spheres_to_add.empty()) {
        return;
    }

    for (const auto &sphere : spheres_to_add) {

        int cell_id = sphere.id;

        if (cell_id < 0) {
            cerr << "Error: Sphere with invalid id: " << cell_id << endl;
            assert(0);
        }

        spheres.push_back(sphere);
        
    }

    skeleton.clear();
    for (unsigned i = 0; i < spheres.size(); i++){
        skeleton.push_back(spheres[i].P);
    } 

    double grid_cell_size = 5e-3;
    double pad = barrier_tickness;
    build_bv_grid_spheres(spheres,grid_cell_size, pad);

}


void Blood_Vessel::build_bv_grid_spheres(const std::vector<Sphere>& spheres_to_add,
                                      double cell_size, double pad)
{

    spheres = spheres_to_add;
    // init grid
    grid = HashGrid{};
    grid.cell = (cell_size > 0.0 ? cell_size : 1.0);
    grid.build_pad = std::max(0.0, pad);

    double maxR = 0.0;
    for (const auto& s : spheres){
        if (s.radius > 0.0){
            maxR = std::max(maxR, s.radius + grid.build_pad);
        }
    }

    grid.max_radius_plus_pad = maxR;

    // 1) compute big box from spheres, already padded by `pad`
    Box B = make_empty_box();

    for (const auto& s : spheres) {
        if (s.radius <= 0.0) continue;
        const double R = s.radius + pad;
        extend(B, s.P - Eigen::Vector3d::Constant(R));
        extend(B, s.P + Eigen::Vector3d::Constant(R));
        // track global max radius+pad
        if (R > grid.max_radius_plus_pad){ 
            grid.max_radius_plus_pad = R;
        }
    }
    

    // 2) set grid origin and store big box
    grid.origin  = Eigen::Vector3d(B.x_min, B.y_min, B.z_min);
    grid.big_box = B;

    // (Optional) reserve to avoid reallocation
    grid.objs.reserve(spheres.size());

    // 3) add spheres into buckets — store indices (b,i)
    auto add = [&](int i) {
        const Sphere& s = spheres[i];
        if (s.radius <= 0.0) return;
        const double R = s.radius + pad;
        const Eigen::Vector3d mn = s.P - Eigen::Vector3d::Constant(R);
        const Eigen::Vector3d mx = s.P + Eigen::Vector3d::Constant(R);

        const Eigen::Array3i imin = ((mn - grid.origin).array() / grid.cell).floor().cast<int>();
        const Eigen::Array3i imax = ((mx - grid.origin).array() / grid.cell).floor().cast<int>();

        const int idx = static_cast<int>(grid.objs.size());
        grid.objs.push_back(i);     // store reference into `spheres`

        for (int ix = imin.x(); ix <= imax.x(); ++ix)
          for (int iy = imin.y(); iy <= imax.y(); ++iy)
            for (int iz = imin.z(); iz <= imax.z(); ++iz)
              grid.buckets[hash3(ix,iy,iz)].push_back(idx);
    };

    for (int i = 0; i < (int)spheres.size(); ++i) {
        const Sphere& s = spheres[i];
        if (s.radius <= 0.0) continue;
        add(i);

    }

}

inline int Blood_Vessel::neighbor_radius_cells(const HashGrid& G, double query_pad)
{
    // We need to cover centers as far as (max sphere radius + query_pad)
    const double max_sphere_radius = std::max(0.0, G.max_radius_plus_pad - G.build_pad);
    const double Rcover = max_sphere_radius + std::max(0.0, query_pad);
    return std::max(1, (int)std::ceil(Rcover / G.cell));
}

inline bool Blood_Vessel::point_in_inflated_aabb(const Eigen::Vector3d& p,
                                   double d)
{
    Box box = grid.big_box;
    const double infl = std::max(0.0, d);
    return (p.x() >= box.x_min - infl && p.x() <= box.x_max + infl &&
            p.y() >= box.y_min - infl && p.y() <= box.y_max + infl &&
            p.z() >= box.z_min - infl && p.z() <= box.z_max + infl);
}


// Optional: segment vs AABB to clamp traversal to the grid big_box (slab method)
inline bool Blood_Vessel::segment_aabb_intersect(const Eigen::Vector3d& p0,
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


void Blood_Vessel::gather_candidates_DDA(const Eigen::Vector3d& p0,
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
inline bool Blood_Vessel::raySphere(const Eigen::Vector3d& p0,
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

bool Blood_Vessel::checkCollision(const Walker& walker,
                           Eigen::Vector3d& step,
                           const double& step_length,
                           Collision& collision)
{

    
    // Normalize direction
    const double L = step_length;

    //cout << "--------- checkCollision Blood_Vessel ---------" << endl;
    //cout <<"L = " << L << endl;

    const Eigen::Vector3d dir = step.normalized();
    const Eigen::Vector3d p0  = walker.pos_v;
    const double Rpad = grid.build_pad;
    bool isinside = isPosInsideBlood_Vessel(p0, barrier_tickness, L);

    if (L <= 0.0) { 
        collision.t = 0.0;
        collision.type = Collision::hit;
        collision.collision_point = walker.pos_v;
        collision.bounced_direction = step.normalized();
        collision.obstacle_type = 3;                    // your code’s convention
        collision.obstacle_ind  = -1;           // or branch/local as needed
        collision.col_location  = isinside ? Collision::inside : Collision::outside;
        collision.perm_crossing = 0.0;

        cout << "problem step length <=0 " << endl;
        //assert(0);
        return true; 
    }


    if (!isinside && walker.location == Walker::intra){

        cout << "problem walker location intra but outside bv " << endl;
        assert(0);

    }


    
    const bool start_inside = (walker.location == Walker::intra);

    // Gather candidates
    std::vector<int> cand_ids;
    gather_candidates_DDA(p0, dir, L, cand_ids);

    struct Ev { double t; int delta; const Sphere* s; };
    std::vector<Ev> evs; evs.reserve(cand_ids.size()*2 + 2);

    const double epsT = std::max(1e-12, 1e-6 * grid.cell);
    
    bool isbouncing = (walker.status == Walker::bouncing);


    auto addSphereEvents = [&](const Sphere* s,
                            const Eigen::Vector3d& p0,
                            const Eigen::Vector3d& dir) {
        if (!s) return;                         
        if (start_inside){

            double t0, t1;
            const double Rin = s->radius;
            if (!raySphere(p0, dir, s->P, Rin, t0, t1)) return;


            // Ensure t0 <= t1 (if your raySphere doesn’t guarantee it)
            if (t1 < t0) std::swap(t0, t1);


            const bool inside0 = (p0 - s->P).squaredNorm() <= (s->radius - Rpad)*(s->radius - Rpad) + 1e-12;
            
            if (!inside0) {
                evs.push_back({ t0, +1, s });
                evs.push_back({ t1, -1, s });
            } 
            else {
                evs.push_back({ t1, -1, s });
            }
        }
        else{

            double t0, t1;
            const double Rin = s->radius;            // ← same inflation here
            if (!raySphere(p0, dir, s->P, Rin, t0, t1)) return;
            if (t1 < t0) std::swap(t0, t1);
            if (t0 > L + Rpad) return;  // intersection beyond step end

            if (t0 >= 0) evs.push_back({t0, +1, s});  // ENTER
            
        }

    };


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
        auto i = grid.objs[idx];   
        const Sphere* sp = &spheres[i];
        if (!sp) { ++null_ptr; std::cerr << "NULL grid.objs["<<idx<<"]\n"; continue; }
        addSphereEvents(sp, p0, dir);
    }
    std::sort(evs.begin(), evs.end(), [](const Ev& a, const Ev& b){ return a.t < b.t; });



    if (bad_idx || null_ptr) {
        std::cerr << "Summary: bad_idx=" << bad_idx << " null_ptr=" << null_ptr << "\n";
    }

    if (evs.empty()) { 
        collision.type = Collision::null; 
        return false; 
    }

    // Drop any events beyond L (in case raySphere or FP noise sneaks one in)
    while (!evs.empty() && evs.back().t > L + 1e-12) evs.pop_back();
    if (evs.empty()) { 

        collision.type = Collision::null; 
        return false; 
    }


    // 5) Sweep to find first union boundary:
    int occ0;

    if (start_inside){
        occ0 = occupancy_at_point(p0, Rpad, start_inside, L);
    }
    else{
        occ0 = occupancy_at_point(p0, -Rpad, start_inside, L);
    }

    bool is_inside = (occ0 > 0);

    int occ = occ0;
    if (start_inside && !is_inside){

        // problem
        collision.type = Collision::hit;
        collision.col_location  = Collision::outside;
        collision.perm_crossing = 0.0;

        //assert(0);
        return true;
    }
    
    else if (!start_inside && is_inside){
        // problem
        collision.type = Collision::hit;
        collision.col_location  = Collision::inside;
        collision.perm_crossing = 0.0;
        //assert(0);
        return true;
    }
    
    
    const Ev* hit = nullptr;

    if (!start_inside) {
        // robust outside path: first time occ becomes > 0
        // (occ was computed earlier; should be 0 here)
        hit  =&evs[0];
    } else {
        // inside -> first time occ becomes 0, grouping same-t events
        size_t i = 0;
        while (i < evs.size()) {
            const double t = evs[i].t;
            int sum = 0;
            const Ev* exitE = nullptr;           // remember any -1 at this t
            size_t j = i;

            // group events with same time (within tolerance)
            while (j < evs.size() && std::fabs(evs[j].t - t) <= epsT) {
                sum += evs[j].delta;
                if (evs[j].delta == -1 && exitE == nullptr) exitE = &evs[j];
                ++j;
            }

            if (occ + sum <= 0) {
                // leaving the union at time t
                hit = exitE ? exitE : &evs[i];   // fall back if all +1 (rare)
                occ += sum;                      // (optional) for logging
                break;
            }

            occ += sum;
            i = j;
        }
    }

    if (!hit) { 

        //cout <<" No hit found\n";
        collision.type = Collision::null; 
        return false; 
    }

    // 6) Fill collision info
    double t_hit = hit->t;

    const Eigen::Vector3d pos = p0 + t_hit * dir;
    
    const Eigen::Vector3d n   = (pos - hit->s->P).normalized();
    const double dn = dir.dot(n);
    const Eigen::Vector3d bounced = dir - 2.0 * dn * n;

    collision.type = Collision::hit;
    collision.collision_point = pos;
    collision.t = t_hit;
    collision.bounced_direction = bounced.normalized();
    collision.obstacle_type = 3;                    // your code’s convention
    collision.obstacle_ind  = hit->s->id;           // or branch/local as needed
    collision.col_location  = start_inside ? Collision::inside : Collision::outside;
    collision.perm_crossing = 0.0;

         
    return true;
}


int Blood_Vessel::occupancy_at_point(const Eigen::Vector3d& p,
                              double margin,
                              const bool& isintra, const double& L) const
{
    // Tiny, scale-aware tiebreaker to avoid boundary chatter:
    // if we *expect* to be inside, be lenient (inflate a hair);
    // if we expect outside, be strict (shrink a hair).

    const double m   = margin;

    int occ = 0;


    // 2) Big-box quick reject for processes (inflate-only; never shrink the box)
    const double infl = std::max(0.0, m);
    const Box& B = grid.big_box;
    if (p.x() < B.x_min - infl || p.x() > B.x_max + infl ||
        p.y() < B.y_min - infl || p.y() > B.y_max + infl ||
        p.z() < B.z_min - infl || p.z() > B.z_max + infl) {
        return occ;
    }

    // 3) Robust neighbor span in grid cells
    double cell_size = grid.cell;
    int number_cells = int(L/cell_size);
    const int Lc = number_cells + 1;

    // 4) Cell index of p
    const Eigen::Array3i ic =
        ((p - grid.origin).array() / grid.cell).floor().cast<int>();

    // 5) Probe neighboring buckets
    std::unordered_set<int> seen; seen.reserve(64);
    for (int dx = -Lc; dx <= Lc; ++dx)
      for (int dy = -Lc; dy <= Lc; ++dy)
        for (int dz = -Lc; dz <= Lc; ++dz) {
            auto it = grid.buckets.find(hash3(ic[0]+dx, ic[1]+dy, ic[2]+dz));
            if (it == grid.buckets.end()) continue;

            for (int idx : it->second) {
                if (!seen.insert(idx).second) continue;

                auto i = grid.objs[idx];
                if ((size_t)i >= spheres.size()) continue;

                const Sphere& s = spheres[i];

                // Allow shrink/inflate by m; skip if degenerate
                const double R = s.radius + m;
                if (R <= 0.0) continue;

                if ((p - s.P).squaredNorm() <= R*R) ++occ;
            }
        }

    return occ;
}


double Blood_Vessel::minDistance(const Walker& w) 
{

    Eigen::Vector3d p = w.pos_v;

    return this->minDistance(p);
}

double Blood_Vessel::minDistance(const Eigen::Vector3d& p){

    if (spheres.empty()) return std::numeric_limits<double>::infinity();
    const Box& box = grid.big_box;

    if (p.x() >= box.x_min && p.x() <= box.x_max &&
        p.y() >= box.y_min && p.y() <= box.y_max &&
        p.z() >= box.z_min && p.z() <= box.z_max) {
        return 0;
    } 

    double dist_x = min(std::abs(p.x() - box.x_min), std::abs(p.x() - box.x_max));
    double dist_y = min(std::abs(p.y() - box.y_min), std::abs(p.y() - box.y_max));
    double dist_z = min(std::abs(p.z() - box.z_min), std::abs(p.z() - box.z_max));
    double minimum = min(dist_x, min(dist_y, dist_z));
    return minimum;

}


bool Blood_Vessel::isPosInsideBlood_Vessel(const Eigen::Vector3d& p, double margin, const double& L) 
{
    const double pad = std::max(0.0, margin);

    const Box& B = grid.big_box;
    if (p.x() < B.x_min - pad || p.x() > B.x_max + pad ||
        p.y() < B.y_min - pad || p.y() > B.y_max + pad ||
        p.z() < B.z_min - pad || p.z() > B.z_max + pad) return false;

    double cell_size = grid.cell;
    int number_cells = int(L/cell_size);
    const int Lc = number_cells + 1;
    const Eigen::Array3i ic = ((p - grid.origin).array() / grid.cell).floor().cast<int>();

    std::unordered_set<int> seen; seen.reserve(64);
    for (int dx=-Lc; dx<=Lc; ++dx)
      for (int dy=-Lc; dy<=Lc; ++dy)
        for (int dz=-Lc; dz<=Lc; ++dz) {
            auto it = grid.buckets.find(hash3(ic[0]+dx, ic[1]+dy, ic[2]+dz));
            if (it == grid.buckets.end()) continue;
            for (int idx : it->second) {
                if (!seen.insert(idx).second) continue;
                auto i = grid.objs[idx];
                const Sphere& s = spheres[i];
                const double R = s.radius + pad;
                if ((p - s.P).squaredNorm() <= R*R + 1e-12) return true;
            }
        }
    return false;
}