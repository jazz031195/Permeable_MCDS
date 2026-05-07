
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
    double min_distance = barrier_tickness;
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

void Blood_Vessel::distance_to_skeleton(const Walker &w, double& min_dist, Eigen::Vector3d& tangent){
    Eigen::Vector3d O = w.pos_v;
    if (skeleton.empty()) {
        min_dist = 0.0; 
        tangent = Eigen::Vector3d::Zero();
        return;
    }
    if (skeleton.size() == 1) {
        min_dist =  (O - skeleton[0]).norm();
        tangent = (O - skeleton[0]).normalized();
        return;
    }
    tangent = (skeleton[skeleton.size()-1] - skeleton[0]).normalized();
    if (w.normal[2] != 0){ 
        tangent = -tangent;
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
    build_bv_grid_spheres(spheres, 5e-3, barrier_tickness);
}

