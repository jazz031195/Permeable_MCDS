#include "obstacle.h"
#include "sphere.h"
#include "constants.h"
#include <math.h>
#include <chrono>

Obstacle::Obstacle():percolation(0.0), prob_cross_e_i(0.0), prob_cross_i_e(0.0), diffusivity_i(0.0), diffusivity_e(0.0), count_perc_crossings(0)
{}

Obstacle::Obstacle(const Obstacle &obs)
{
    id = obs.id;
    spheres = obs.spheres;
    grid = obs.grid;
    percolation = obs.percolation;
    prob_cross_e_i = obs.prob_cross_e_i;
    prob_cross_i_e = obs.prob_cross_i_e;
    diffusivity_i = obs.diffusivity_i;
    diffusivity_e = obs.diffusivity_e;
    count_perc_crossings = obs.count_perc_crossings;
}

void Obstacle::setPercolation(double &percolation_)
{
    percolation = percolation_;
}

void Obstacle::setDiffusion(double &diffusivity_i_, double &diffusivity_e_){
    diffusivity_i = diffusivity_i_;
    diffusivity_e = diffusivity_e_;
}


void Obstacle::setProbabilities(double &prob_cross_e_i_, double &prob_cross_i_e_)
{
    prob_cross_e_i = prob_cross_e_i_;
    prob_cross_i_e = prob_cross_i_e_;
}


void Obstacle::set_prob_crossings(double step_length_pref){

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
    for (unsigned j = 0; j < spheres.size(); j++){
        spheres[j].percolation = this->percolation;
        if(spheres[j].percolation > 0.0){

            spheres[j].prob_cross_e_i = this->prob_cross_e_i;
            spheres[j].prob_cross_i_e = this->prob_cross_i_e;

            spheres[j].diffusivity_e = this->diffusivity_e;
            spheres[j].diffusivity_i = this->diffusivity_i;
        }
    }
}

void Obstacle::elasticBounceAgainsPlane(Eigen::Vector3d &ray_origin, Eigen::Vector3d &normal, double &t, Eigen::Vector3d &step)
{

    Eigen::Vector3d ray =  (-t*step).normalized();//
    double rn = ray.dot(normal);

    //std::cout << "rn in elastic plane :" << rn << std::endl;

    // Caso 3) ni cerca ni paralela
    step = -ray + 2.0*normal*rn;

    //step = (rn>0.0)?normal:(-normal);

}

double Obstacle::clamp(double value, double lower, double upper) {
    return std::max(lower, std::min(value, upper));
}

inline bool Obstacle::is_empty(const Box& b) {
    return b.x_min > b.x_max || b.y_min > b.y_max || b.z_min > b.z_max;
}

inline void Obstacle::extend(Box& b, const Eigen::Vector3d& p) {
    if (is_empty(b)) { b = {p.x(),p.x(),p.y(),p.y(),p.z(),p.z()}; return; }
    b.x_min = std::min(b.x_min, p.x()); b.x_max = std::max(b.x_max, p.x());
    b.y_min = std::min(b.y_min, p.y()); b.y_max = std::max(b.y_max, p.y());
    b.z_min = std::min(b.z_min, p.z()); b.z_max = std::max(b.z_max, p.z());
}

inline uint64_t hash3(int x, int y, int z) {
    auto mix = [](uint64_t v){ v += 0x9e3779b97f4a7c15ULL; v = (v^(v>>30))*0xbf58476d1ce4e5b9ULL;
                               v = (v^(v>>27))*0x94d049bb133111ebULL; return v^(v>>31); };
    return mix((uint64_t)(uint32_t)x) ^ (mix((uint64_t)(uint32_t)y)<<1) ^ (mix((uint64_t)(uint32_t)z)<<2);
}

// =========================================================================
// GRID BUILDER (SIMPLIFIED)
// =========================================================================


void Obstacle::set_spheres(std::vector<Sphere> &spheres_to_add) {
    spheres.clear();
    if (spheres_to_add.empty()) return;

    for (const auto &sphere : spheres_to_add) {
        if (sphere.id < 0) assert(0);
        spheres.push_back(sphere);
    }
    build_bv_grid_spheres(spheres, 1e-2, barrier_tickness);
}


void Obstacle::build_bv_grid_spheres(const std::vector<Sphere>& spheres_to_add, double cell_size, double pad)
{
    spheres = spheres_to_add;
    grid = HashGrid{};
    grid.cell = (cell_size > 0.0 ? cell_size : 1.0);
    grid.build_pad = std::max(0.0, pad);

    double maxR = 0.0;
    Box B = make_empty_box();

    for (const auto& s : spheres) {
        if (s.radius <= 0.0) continue;
        const double R = s.radius + pad;
        extend(B, s.P - Eigen::Vector3d::Constant(R));
        extend(B, s.P + Eigen::Vector3d::Constant(R));
        maxR = std::max(maxR, R);
    }
    
    grid.max_radius_plus_pad = maxR;
    grid.origin  = Eigen::Vector3d(B.x_min, B.y_min, B.z_min);
    grid.big_box = B;
    grid.objs.reserve(spheres.size());

    for (int i = 0; i < (int)spheres.size(); ++i) {
        const Sphere& s = spheres[i];
        if (s.radius <= 0.0) continue;
        
        const double R = s.radius + pad;
        const Eigen::Vector3d mn = s.P - Eigen::Vector3d::Constant(R);
        const Eigen::Vector3d mx = s.P + Eigen::Vector3d::Constant(R);

        const Eigen::Array3i imin = ((mn - grid.origin).array() / grid.cell).floor().cast<int>();
        const Eigen::Array3i imax = ((mx - grid.origin).array() / grid.cell).floor().cast<int>();

        int idx = static_cast<int>(grid.objs.size());
        grid.objs.push_back(i);     

        for (int ix = imin.x(); ix <= imax.x(); ++ix)
          for (int iy = imin.y(); iy <= imax.y(); ++iy)
            for (int iz = imin.z(); iz <= imax.z(); ++iz)
              grid.buckets[hash3(ix,iy,iz)].push_back(idx);
    }
}


// =========================================================================
// CLEAN COLLISION & INTERSECTION LOGIC
// =========================================================================

// Simple math intersection
bool Obstacle::raySphere(const Eigen::Vector3d& p0, const Eigen::Vector3d& dir_unit, const Eigen::Vector3d& C, double R, double& t_enter, double& t_exit)
{
    const Eigen::Vector3d oc = p0 - C;
    const double b = oc.dot(dir_unit);
    const double c = oc.squaredNorm() - R*R;
    const double disc = b*b - c;
    if (disc < 0.0) return false;

    const double s = std::sqrt(disc);
    t_enter = -b - s;    
    t_exit  = -b + s;    
    if (t_enter > t_exit) std::swap(t_enter, t_exit);
    return true;
}
bool Obstacle::isPosInsideObstacle(const Eigen::Vector3d& p, double margin) 
{
    // ========================================================================
    // GUARD 0: The NaN Shield 
    // Prevents casting NaN to int, which causes Undefined Behavior and segfaults
    // ========================================================================
    if (std::isnan(p.x()) || std::isnan(p.y()) || std::isnan(p.z())) {
        std::cout << "\n[FATAL ERROR] NaN position passed to isPosInsideObstacle!" << std::endl;
        assert(0);
        return false;
    }

    // 1. Fast-Fail against the Global Bounding Box
    // It is mathematically safe to keep 'margin' here as a search buffer.
    const Box& B = grid.big_box;
    if (p.x() < B.x_min - margin || p.x() > B.x_max + margin ||
        p.y() < B.y_min - margin || p.y() > B.y_max + margin ||
        p.z() < B.z_min - margin || p.z() > B.z_max + margin) {
        return false;
    }

    // 2. Map position to grid buckets safely
    const double inv_cell = 1.0 / grid.cell; 
    const int Lc = static_cast<int>(grid.max_radius_plus_pad * inv_cell) + 1;
    const Eigen::Array3i ic = ((p - grid.origin).array() * inv_cell).floor().cast<int>();

    // 3. Search neighboring buckets
    for (int dx = -Lc; dx <= Lc; ++dx) {
        for (int dy = -Lc; dy <= Lc; ++dy) {
            for (int dz = -Lc; dz <= Lc; ++dz) {
                
                auto it = grid.buckets.find(hash3(ic[0] + dx, ic[1] + dy, ic[2] + dz));
                if (it == grid.buckets.end()) continue;

                // 4. Safely evaluate candidates
                for (int idx : it->second) {
                    
                    // ========================================================
                    // GUARD 1: Prevent grid.objs out-of-bounds access
                    // ========================================================
                    if (idx < 0 || idx >= grid.objs.size()) {
                        std::cout << "\n[FATAL MEMORY CORRUPTION] grid.objs size is " << grid.objs.size() 
                                  << ", but spatial hash requested idx: " << idx << std::endl;
                        assert(0);
                    }

                    int sphere_id = grid.objs[idx];

                    // ========================================================
                    // GUARD 2: Prevent spheres array out-of-bounds access
                    // ========================================================
                    if (sphere_id < 0 || sphere_id >= spheres.size()) {
                        std::cout << "\n[FATAL MEMORY CORRUPTION] spheres vector size is " << spheres.size() 
                                  << ", but grid points to sphere_id: " << sphere_id << std::endl;
                        assert(0);
                    }

                    // ========================================================
                    // 5. THE FIX: Strict Physical Intersection Math
                    // ========================================================
                    const Sphere& s = spheres[sphere_id];
                    
                    // We DO NOT add margin to the radius anymore.
                    // We subtract a tiny epsilon to protect against floating-point rounding errors
                    // incorrectly identifying edge-walkers as "inside".
                    const double epsilon = 1e-12; 
                    const double true_squared_radius = (s.radius * s.radius) + margin;
                    
                    if ((p - s.P).squaredNorm() <= true_squared_radius) {
                        return true; // Point is strictly inside this sphere!
                    }
                }
            }
        }
    }
    
    // Checked all candidates in range, no intersection found
    return false;
}
void Obstacle::gather_candidates_AABB(const Eigen::Vector3d& p0, const Eigen::Vector3d& dir, double L, const HashGrid& grid, std::vector<int>& out_ids)
{
    out_ids.clear();
    Eigen::Vector3d p1 = p0 + dir * L;
    
    // Strict bounding box of the ray. No inflation needed!
    Eigen::Vector3d min_pt = p0.cwiseMin(p1);
    Eigen::Vector3d max_pt = p0.cwiseMax(p1);

    const double inv_cell = 1.0 / grid.cell;
    Eigen::Array3i imin = ((min_pt - grid.origin).array() * inv_cell).floor().cast<int>();
    Eigen::Array3i imax = ((max_pt - grid.origin).array() * inv_cell).floor().cast<int>();

    for(int ix = imin.x(); ix <= imax.x(); ++ix) {
        for(int iy = imin.y(); iy <= imax.y(); ++iy) {
            for(int iz = imin.z(); iz <= imax.z(); ++iz) {
                auto it = grid.buckets.find(hash3(ix,iy,iz));
                if (it != grid.buckets.end()) {
                    for (int idx : it->second) {
                        out_ids.push_back(idx);
                    }
                }
            }
        }
    }

    // Fast deduplication
    if (!out_ids.empty()) {
        std::sort(out_ids.begin(), out_ids.end());
        out_ids.erase(std::unique(out_ids.begin(), out_ids.end()), out_ids.end());
    }
}

bool Obstacle::checkCollision(const Walker& walker, Eigen::Vector3d& step, const double& step_length, Collision& collision)
{

    const double L = step_length;
    const Eigen::Vector3d p0  = walker.pos_v;

    if (L <= 0.0) { 
        collision.type = Collision::null; 
        cout << "\n[FATAL ERROR] Step length must be positive!" << endl;
        assert(0); // Should never happen, step length must be positive
        return false; 
    }

    // check if p0 is nan
    if (std::isnan(p0.x()) || std::isnan(p0.y()) || std::isnan(p0.z())) {
        cout << "\n[FATAL ERROR] Walker position is NaN!" << endl;
        assert(0);
    }

    if (step.squaredNorm() < 1e-14) {
        cout << "\n[FATAL ERROR] step is smaller than 1e-14 !" << endl;
        cout << "step: " << step.transpose() << endl;
        assert(0);
    }

    if (std::isnan(step.x()) || std::isnan(step.y()) || std::isnan(step.z())) {
        cout << "\n[FATAL ERROR] step is NaN!" << endl;
        assert(0);
    }

    const Eigen::Vector3d dir = step.normalized();

    // 1. Establish absolute truth
    bool math_inside_strict = isPosInsideObstacle(p0, barrier_tickness);
    //bool math_inside_strict = true;
    bool start_inside = (walker.location == Walker::intra);

    // ====================================================================
    // STRICT LEAK CHECKER (The Tolerance Trap)
    // ====================================================================
    if (math_inside_strict != start_inside) {
   
        // 1e-5 mm tolerance strictly for forgiving harmless 64-bit float drift 
        // that occurs naturally when resting directly on a sphere boundary.
        double drift_tolerance = 1e-5; 

        if (start_inside && !math_inside_strict) {
            if (isPosInsideObstacle(p0, barrier_tickness + drift_tolerance)) {
                start_inside = true; // Forgive micro-drift, trust the walker
            } else {
                cout << "\n[FATAL LEAK] Walker state: 'intra', but physically OUTSIDE the vessel!" << endl;
                cout << "Walker Pos: " << p0.transpose() << endl;
                collision.type = Collision::leaked; 
                return false; 
            }
        }
        else if (!start_inside && math_inside_strict) {
            if (!isPosInsideObstacle(p0, barrier_tickness - drift_tolerance)) {
                start_inside = false; // Forgive micro-drift, trust the walker
            } else {
                cout << "\n[FATAL LEAK] Walker state: 'extra', but physically INSIDE the vessel!" << endl;
                cout << "Walker Pos: " << p0.transpose() << endl;
                collision.type = Collision::leaked; 
                return false; 
            }
        }
    }

    // 2. Gather Candidates based on ray bounding box
    static thread_local std::vector<int> cand_ids;
    gather_candidates_AABB(p0, dir, L, grid, cand_ids);

    if (cand_ids.empty()) {
        collision.type = Collision::null;

        return false;
    }

    // 3. Find all potential intersections
    struct Hit { double t; const Sphere* s; bool is_exit; };
    static thread_local std::vector<Hit> hits;
    hits.clear(); // Must clear it manually before use!
    hits.reserve(cand_ids.size() * 2);


    for (int idx : cand_ids) {
        const Sphere& s = spheres[grid.objs[idx]];
        double t0, t1;
        
        if (raySphere(p0, dir, s.P, s.radius, t0, t1)) {
            // No blind spots. Keep all mathematical hits ahead of the ray up to L.
            if (t0 >= -1e-9 && t0 <= L + 1e-9) hits.push_back({std::max(0.0, t0), &s, false});
            if (t1 >= -1e-9 && t1 <= L + 1e-9) hits.push_back({std::max(0.0, t1), &s, true});
        }
    }


    if (hits.empty()) {
        collision.type = Collision::null;

        return false;
    }

    // Sort hits chronologically
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b){ return a.t < b.t; });

    // 4. CSG Boundary Evaluation
    const Hit* valid_hit = nullptr;

    for (const auto& hit : hits) {
        Eigen::Vector3d pt = p0 + hit.t * dir;
        bool buried = false;

        // Check if this intersection point is physically "inside" any OTHER sphere
        for (int idx : cand_ids) {
            const Sphere& other_s = spheres[grid.objs[idx]];
            if (&other_s == hit.s) continue;
            
            // Shrink other sphere slightly to ignore surface seams cleanly
            const double shrunk_radius = std::max(0.0, other_s.radius - 1e-7);
            if ((pt - other_s.P).squaredNorm() <= (shrunk_radius * shrunk_radius)) {
                buried = true;
                break; // Buried inside the union! Ignore this boundary.
            }
        }

        if (!buried) {
            // The true, unburied boundary of the Union
            if (start_inside && hit.is_exit) { 
                valid_hit = &hit; 
                break; 
            }
            if (!start_inside && !hit.is_exit) { 
                valid_hit = &hit; 
                break; 
            }
        }
    }
    

    if (!valid_hit) { 
        collision.type = Collision::null; 

        return false; 
    }

    // 5. Build final collision response
    collision.type = Collision::hit;
 
    if (std::isnan(valid_hit->t)) {
        cout << "\n[FATAL ERROR] valid_hit->t is NaN!" << endl;
        assert(0);
    }
    collision.t = valid_hit->t;
    collision.collision_point = p0 + valid_hit->t * dir;
    
    Eigen::Vector3d n = (collision.collision_point - valid_hit->s->P).normalized();
    const double dn = dir.dot(n);
    collision.bounced_direction = (dir - 2.0 * dn * n).normalized();
    
    collision.obstacle_type = this->getObstacleType();              
    collision.obstacle_ind  = valid_hit->s->id;           
    collision.col_location  = start_inside ? Collision::inside : Collision::outside;
    collision.perm_crossing = 0.0;

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
// Misc utility (Untouched)
double Obstacle::minDistance(const Walker& w) { return this->minDistance(w.pos_v); }

double Obstacle::minDistance(const Eigen::Vector3d& p){
    if (spheres.empty()) return std::numeric_limits<double>::infinity();
    const Box& box = grid.big_box;
    if (p.x() >= box.x_min && p.x() <= box.x_max && p.y() >= box.y_min && p.y() <= box.y_max && p.z() >= box.z_min && p.z() <= box.z_max) {
        return 0;
    } 
    double dist_x = min(std::abs(p.x() - box.x_min), std::abs(p.x() - box.x_max));
    double dist_y = min(std::abs(p.y() - box.y_min), std::abs(p.y() - box.y_max));
    double dist_z = min(std::abs(p.z() - box.z_min), std::abs(p.z() - box.z_max));
    return min(dist_x, min(dist_y, dist_z));
}
