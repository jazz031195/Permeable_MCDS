
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
#include <algorithm> // for std::min, std::unique, std::sort

using namespace Eigen;
using namespace std;

std::mt19937 gen;

Blood_Vessel::Blood_Vessel() {}

Blood_Vessel::~Blood_Vessel() {}

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
    min_velocity = bv.min_velocity;
    use_blood_random_direction = bv.use_blood_random_direction;
};

int Blood_Vessel::getObstacleType() const { 
    return blood_obstacle_type; 
}


// =========================================================================
// FLOW, VELOCITY, & SKELETON (UNTOUCHED)
// =========================================================================

void Blood_Vessel::set_bv_parameters(const double &pressure_diff_){
    pressure_diff = pressure_diff_;
    flow = (M_PI*pressure_diff*(radius*radius*radius*radius))/(8*viscosity); 
    max_velocity = (pressure_diff/(4*viscosity))*(radius*radius); 
    //cout << "radius :" << radius << ", max velocity: " << max_velocity << endl;
    double min_distance = 1e-6;
    min_velocity = velocity(radius - min_distance);
}

double distancePointToSegment(const Eigen::Vector3d& P, const Eigen::Vector3d& A, const Eigen::Vector3d& B) {
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AP = P - A;
    double ab2 = (AB).norm() * (AB).norm();
    if (ab2 == 0.0) return (P - A).norm();
    double t = (AP).dot(AB) / ab2;
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;
    Eigen::Vector3d closest = A + AB * t;
    return (P - closest).norm();
}
void Blood_Vessel::distance_to_skeleton(const Walker &w, double& min_dist, Eigen::Vector3d& tangent) {
    Eigen::Vector3d O = w.pos_v;

    if (std::isnan(O[0])) {
        std::cout << "CRASH: Walker " << w.index << " is NaN BEFORE distance_to_skeleton!" << std::endl;
        assert(0); 
    }

    if (skeleton.empty() || skeleton.size() == 1) {
        min_dist = 0.0; 
        tangent = Eigen::Vector3d::Zero();
        assert(0 && "Blood_Vessel skeleton is invalid!");
        return;
    }

    double min_sq_dist = std::numeric_limits<double>::infinity();
    Eigen::Vector3d best_tangent = Eigen::Vector3d::Zero();

    // Loop through every connected segment
    for (size_t i = 0; i < skeleton.size() - 1; ++i) {
        const Eigen::Vector3d& A = skeleton[i];
        const Eigen::Vector3d& B = skeleton[i + 1];

        Eigen::Vector3d AB = B - A;
        double length_sq = AB.squaredNorm();

        if (length_sq < 1e-12) continue;

        double t = (O - A).dot(AB) / length_sq;

        // ==========================================
        // FIX 1: THE ANTI-SHRINK BOUNDARY EXTRAPOLATION
        // ==========================================
        if (i == 0) {
            t = std::min(1.0, t); // First segment: Extend infinitely backward
        } else if (i == skeleton.size() - 2) {
            t = std::max(0.0, t); // Last segment: Extend infinitely forward
        } else {
            t = std::max(0.0, std::min(1.0, t)); // Internal segments: Clamp strictly
        }

        Eigen::Vector3d closest_pt_on_segment = A + t * AB;
        double sq_dist = (O - closest_pt_on_segment).squaredNorm();

        if (sq_dist < min_sq_dist) {
            min_sq_dist = sq_dist;
            best_tangent = AB.normalized();
        }
    }

    if (std::isinf(min_sq_dist)) {
        min_dist = (O - skeleton.front()).norm();
        tangent = Eigen::Vector3d::Zero();
        assert(0);
    } else {
        min_dist = std::sqrt(min_sq_dist);
        tangent = best_tangent;
    }

    // ==========================================
    // FIX 2: THE LOOPING/FLIPPING TRAP
    // ==========================================
    // If you are using Periodic Boundaries (wrapping) for your tortuous vessels,
    // the walker flows continuously. Flipping the tangent here will instantly 
    // force the walker to turn around and flow backward, causing infinite loops!
    // 
    // ONLY uncomment these lines if you are simulating a perfectly straight, 
    // mathematically symmetric pipe with mirror walls.
    
    // 5. Apply the Unified Mirror Flow Toggle (The Parity Bit)
    // Sum the mirror states to determine if we are flowing forward or backward
    int flip_count = w.normal[0] + w.normal[1] + w.normal[2];
    
    // If an odd number of axes are mirrored, we must flow exactly backward 
    // down all 3 dimensions of the physical pipe
    if (flip_count % 2 != 0) {
        tangent = -tangent;
    }

    if (tangent.squaredNorm() > 1e-12) {
        tangent.normalize();
    }
    

    if (tangent.squaredNorm() > 1e-12) {
        tangent.normalize();
    }
}

double Blood_Vessel::velocity(const double &radial_distance){
    return (pressure_diff/(4*viscosity))*(this->radius*this->radius - radial_distance*radial_distance);
}

void Blood_Vessel::WalkerVelocity(const Walker &w, double& v, Eigen::Vector3d& flow_direction){
    double min_dist;
    Eigen::Vector3d tangent;
    distance_to_skeleton(w, min_dist, tangent);

    if (min_dist >= this->radius) {
        v = min_velocity;
        flow_direction = tangent;
        return;
    }
    double radial_distance = min_dist; 
    v = velocity(radial_distance);
    if (v < min_velocity){
        v = min_velocity;
    }
    if (use_blood_random_direction){
        flow_direction = biased_direction_from_tangent(tangent);
    } else {
        flow_direction = tangent;
    }
}

Eigen::Vector3d Blood_Vessel::biased_direction_from_tangent(const Eigen::Vector3d& tangent) {
    Eigen::Vector3d target_direction = tangent.normalized();
    double std_dev = 0.1; 
    Eigen::Vector3d direction = apply_bias_toward_target(
        generate_random_point_on_sphere(std_dev),
        target_direction
    ).normalized();
    return direction;
}

Eigen::Vector3d Blood_Vessel::apply_bias_toward_target(const Eigen::Vector3d &point, const Eigen::Vector3d &target) {
    Eigen::Vector3d reference(0, 0, 1); 
    Eigen::Matrix3d R = rotation_matrix_from_vectors(reference, target);
    Eigen::Vector3d rotated_point = (R * point).normalized();
    return rotated_point;
}

Eigen::Vector3d Blood_Vessel::generate_random_point_on_sphere(double std) {
    if (std == 0) {
        return Eigen::Vector3d(0, 0, 1);
    }
    std::normal_distribution<double> N(0.0, std);
    double phi = std::abs(N(gen));
    phi = std::min(phi, M_PI / 2.0);  
    std::uniform_real_distribution<double> U(0.0, 2.0 * M_PI);
    double theta = U(gen);
    double x = std::sin(phi) * std::cos(theta);
    double y = std::sin(phi) * std::sin(theta);
    double z = std::cos(phi);
    return Eigen::Vector3d(x, y, z);  
}

Eigen::Matrix3d Blood_Vessel::rotation_matrix_from_vectors(const Eigen::Vector3d &vec1, const Eigen::Vector3d &vec2) {
    Eigen::Vector3d a = vec1.normalized();
    Eigen::Vector3d b = vec2.normalized();
    double c = a.dot(b);
    if (c > 1.0 - 1e-12) return Eigen::Matrix3d::Identity();       
    if (c < -1.0 + 1e-12) {                                         
        Eigen::Vector3d axis = a.unitOrthogonal();
        return Eigen::AngleAxisd(M_PI, axis).toRotationMatrix();
    }
    Eigen::Vector3d v = a.cross(b);
    double s = v.norm();
    Eigen::Matrix3d K;
    K <<   0,   -v.z(),  v.y(),
         v.z(),     0,  -v.x(),
        -v.y(),  v.x(),    0;
    return Eigen::Matrix3d::Identity() + K + K*K * ((1 - c)/(s*s));
}


// =========================================================================
// GRID BUILDER (SIMPLIFIED)
// =========================================================================

void Blood_Vessel::set_spheres(std::vector<Sphere> &spheres_to_add) {
    spheres.clear();
    if (spheres_to_add.empty()) return;

    for (const auto &sphere : spheres_to_add) {
        if (sphere.id < 0) assert(0);
        spheres.push_back(sphere);
    }
    skeleton.clear();
    for (unsigned i = 0; i < spheres.size(); i++){
        skeleton.push_back(spheres[i].P);
    } 
    build_bv_grid_spheres(spheres, 1, barrier_tickness);
}

bool Blood_Vessel::checkCollision(const Walker& walker, Eigen::Vector3d& step, const double& step_length, Collision& collision)
{
    const double L = step_length;
    const Eigen::Vector3d p0  = walker.pos_v;

    // 1. Basic Safety Guards
    if (L <= 0.0) { 
        collision.type = Collision::null; 
        return false; 
    }
    if (std::isnan(p0.x()) || std::isnan(p0.y()) || std::isnan(p0.z()) || 
        std::isnan(step.x()) || std::isnan(step.y()) || std::isnan(step.z())) {
        cout << "\n[FATAL ERROR] NaN detected in position or step!" << endl;
        assert(0);
    }
    if (step.squaredNorm() < 1e-14) {
        cout << "\n[FATAL ERROR] step is smaller than 1e-14 !" << endl;
        assert(0);
    }

    const Eigen::Vector3d dir = step.normalized();
    bool start_inside = (walker.location == Walker::intra);

    // ====================================================================
    // THE INTRAVASCULAR FAST-PATH (Massive Simplification)
    // ====================================================================
    // Since blood flow follows the vessel, 99.9% of steps will remain inside.
    // If the destination is still inside the vessel network, skip ALL complex math!
    if (start_inside) {
        Eigen::Vector3d p_end = p0 + step;
        
        if (isPosInsideObstacle(p_end, barrier_tickness)) {
            // Walker flowed safely down the pipe. No boundary hit!
            collision.type = Collision::null;
            return false; 
        }
        
        // If we reach here, it means the walker actually hit the wall 
        // (e.g., due to diffusion or a sharp tortuous curve).
        // We let it fall through to the CSG logic below to bounce or permeate.
    }

    // ====================================================================
    // STRICT LEAK CHECKER (Only runs if we might hit a boundary)
    // ====================================================================
    bool math_inside_strict = isPosInsideObstacle(p0, barrier_tickness);
    
    if (math_inside_strict != start_inside) {
        double drift_tolerance = 1e-5; 
        if (start_inside && !math_inside_strict) {
            if (isPosInsideObstacle(p0, barrier_tickness + drift_tolerance)) {
                start_inside = true; 
            } else {
                collision.type = Collision::leaked; 
                return false; 
            }
        }
        else if (!start_inside && math_inside_strict) {
            if (!isPosInsideObstacle(p0, barrier_tickness - drift_tolerance)) {
                start_inside = false; 
            } else {
                collision.type = Collision::leaked; 
                return false; 
            }
        }
    }

    // ====================================================================
    // STANDARD CSG RAYCASTING (For tissue walkers & blood bouncing)
    // ====================================================================
    static thread_local std::vector<int> cand_ids;
    gather_candidates_AABB(p0, dir, L, grid, cand_ids);

    if (cand_ids.empty()) {
        collision.type = Collision::null;
        return false;
    }

    struct Hit { double t; const Sphere* s; bool is_exit; };
    static thread_local std::vector<Hit> hits;
    hits.clear(); 
    hits.reserve(cand_ids.size() * 2);

    for (int idx : cand_ids) {
        const Sphere& s = spheres[grid.objs[idx]];
        double t0, t1;
        
        if (raySphere(p0, dir, s.P, s.radius, t0, t1)) {
            if (t0 >= -1e-9 && t0 <= L + 1e-9) hits.push_back({std::max(0.0, t0), &s, false});
            if (t1 >= -1e-9 && t1 <= L + 1e-9) hits.push_back({std::max(0.0, t1), &s, true});
        }
    }

    if (hits.empty()) {
        collision.type = Collision::null;
        return false;
    }

    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b){ return a.t < b.t; });

    const Hit* valid_hit = nullptr;

    for (const auto& hit : hits) {
        Eigen::Vector3d pt = p0 + hit.t * dir;
        bool buried = false;

        for (int idx : cand_ids) {
            const Sphere& other_s = spheres[grid.objs[idx]];
            if (&other_s == hit.s) continue;
            
            const double shrunk_radius = std::max(0.0, other_s.radius - 1e-7);
            if ((pt - other_s.P).squaredNorm() <= (shrunk_radius * shrunk_radius)) {
                buried = true;
                break; 
            }
        }

        if (!buried) {
            if (start_inside && hit.is_exit) { valid_hit = &hit; break; }
            if (!start_inside && !hit.is_exit) { valid_hit = &hit; break; }
        }
    }

    if (!valid_hit) { 
        collision.type = Collision::null; 
        return false; 
    }

    // ====================================================================
    // BUILD COLLISION RESPONSE
    // ====================================================================
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

    if (percolation > 0.0) {
        static thread_local std::mt19937 gen{std::random_device{}()};
        std::uniform_real_distribution<double> U(0.0,1.0);
        const double p_cross = start_inside ? prob_cross_i_e : prob_cross_e_i;
        if (U(gen) < p_cross) {
            collision.perm_crossing = p_cross;
            collision.bounced_direction = dir; 
        }
    }

    return true;
}