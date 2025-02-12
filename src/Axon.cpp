#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include <numeric> // Add this line to include std::iota
#include <unordered_set>
#include <iomanip>

using namespace Eigen;
using namespace std;


Axon::Axon(const Axon &ax)
{
    id = ax.id;
    spheres = ax.spheres;
    radius = ax.radius;
    begin = ax.begin;
    end = ax.end;
    boxes = ax.boxes;
    percolation = ax.percolation;
    prob_cross_e_i = ax.prob_cross_e_i;
    prob_cross_i_e = ax.prob_cross_i_e;
    diffusivity_i = ax.diffusivity_i;
    diffusivity_e = ax.diffusivity_e;

}

void Axon::set_spheres(const std::vector<Sphere> &spheres_to_add) {

    boxes.clear();
    spheres.clear();

    if (spheres_to_add.empty()) {
        std::cout << "No spheres to add." << std::endl;
        return;
    }

    int nbr_spheres_per_box = spheres_to_add.size() / nbr_boxes;
    int remainder = spheres_to_add.size() % nbr_boxes;

    boxes.resize(nbr_boxes);  // Pre-construct the required number of boxes

    int index = 0;

    for (int i = 0; i < nbr_boxes; ++i) {
        Box &box = boxes[i];
        box.x_min = box.y_min = box.z_min = std::numeric_limits<double>::max();
        box.x_max = box.y_max = box.z_max = std::numeric_limits<double>::lowest();

        int spheres_in_this_box = nbr_spheres_per_box + (i < remainder ? 1 : 0);

        for (int j = 0; j < spheres_in_this_box && index < spheres_to_add.size(); ++j, ++index) {
            const Sphere &sphere = spheres_to_add[index];

            // Update object ID and sphere ID
            spheres.push_back(sphere);
            spheres.back().object_id = id;
            spheres.back().id = index;

            // Initialize or expand box boundaries
            box.x_min = std::min(box.x_min, sphere.P[0] - sphere.radius);
            box.x_max = std::max(box.x_max, sphere.P[0] + sphere.radius);
            box.y_min = std::min(box.y_min, sphere.P[1] - sphere.radius);
            box.y_max = std::max(box.y_max, sphere.P[1] + sphere.radius);
            box.z_min = std::min(box.z_min, sphere.P[2] - sphere.radius);
            box.z_max = std::max(box.z_max, sphere.P[2] + sphere.radius);
        }
    }

    // Set the beginning and end positions
    this->begin = spheres_to_add.front().P;
    this->end = spheres_to_add.back().P;

    // Compute mean radius
    double mean_radius = 0.0;
    for (const auto &sphere : spheres) {
        mean_radius += sphere.radius;
    }
    this->radius = spheres.empty() ? 0.0 : mean_radius / spheres.size();

}

bool Axon::intersection_sphere_vector(double &t1, double &t2, const Sphere &s, const Eigen::Vector3d &step, const Eigen::Vector3d &pos) {
    Eigen::Vector3d m = pos - s.P;
    double rad2 = s.radius * s.radius;

    double b = m.dot(step);
    double c = m.squaredNorm() - rad2;
    double discr = b * b - c;

    if (discr < 0.0) {
        return false; // No intersection
    }

    double sqrt_discr = std::sqrt(discr);
    double neg_b = -b;  // Avoid recalculating -b twice

    t1 = neg_b + sqrt_discr;
    t2 = neg_b - sqrt_discr;

    return true;
}

void Axon::find_all_intersections(const Walker &walker, const Eigen::Vector3d &step, const double &distance,
                                  std::vector<std::pair<double, size_t>> &dist_and_indices) {
    dist_and_indices.clear();

    std::vector<int> sph_ids_walker_is_inside;
    Eigen::Vector3d pos = walker.pos_v;

    // Find spheres the walker is near
    if (!FindSphereinAxon(pos, distance, sph_ids_walker_is_inside)) {
        //cout << "not near axon" << endl;
        return;
    }

    dist_and_indices.reserve(sph_ids_walker_is_inside.size());

    // Process each sphere the walker is near
    for (int sphere_id : sph_ids_walker_is_inside) {
        Sphere sphere_to_check = spheres[sphere_id];

        double t1, t2;
        if (intersection_sphere_vector(t1, t2, sphere_to_check, step, pos)) {
            // Check both intersections
            auto check_and_add = [&](double t) {
                if (t > ((walker.status == Walker::bouncing) ? EPS_VAL : 0)) {
                    dist_and_indices.emplace_back(t, sphere_id);
                }
            };
            check_and_add(t1);
            check_and_add(t2);
        }
    }
}

bool Axon::checkCollision(const Walker &walker, Eigen::Vector3d &step, const double &step_length, Collision &collision) {
    
    // Distances to intersections and corresponding sphere IDs
    std::vector<std::pair<double, size_t>> dist_and_indices;

    // Find all intersections
    find_all_intersections(walker, step, step_length + barrier_tickness, dist_and_indices);

    if (dist_and_indices.empty()) {
        // Handle case with no intersections
        
        if (walker.location == Walker::intra) {

            if (!isPosInsideAxon_(walker.pos_v, EPS_VAL)){ 

                //cout << "position : " << walker.pos_v << endl;
                //cout << "is not inside id : " << id << endl;
                collision.type = Collision::hit;
                collision.col_location = Collision::outside;
                collision.collision_point = walker.pos_v;
                collision.t = 1e-9;
                collision.perm_crossing = 0.0;
                collision.bounced_direction = step;

                return true;
            } 
        }
        
        collision.type = Collision::null;
        return false;
    }
    if (walker.location == Walker::intra){
        std::sort(dist_and_indices.begin(), dist_and_indices.end());
    }
    else{
        // Find the closest distance (smallest element)
        std::nth_element(dist_and_indices.begin(), dist_and_indices.begin(), dist_and_indices.end());
    }  

    for (auto dist_and_indice : dist_and_indices) {
        
        size_t i = dist_and_indice.second;
        double distance = dist_and_indice.first;

        if (distance > step_length + barrier_tickness) {
            continue; // Early exit: remaining distances are too far
        }

        Eigen::Vector3d pos = walker.pos_v + distance * step;

        // Check if position is near the edge
        bool is_near_edge = (walker.location == Walker::intra) ? !isPosInsideAxon_(pos, -EPS_VAL) : true;

        if (is_near_edge && distance < step_length + barrier_tickness) {

            Sphere sphere = spheres[i];

            // Compute normal and bounced direction
            Eigen::Vector3d normal = (pos - sphere.P).normalized();
            Eigen::Vector3d ray = (-distance * step).normalized();
            double rn = ray.dot(normal);

            if (rn < -1e-10){
                collision.col_location = Collision::inside;
            }
            //outside
            else if (rn > 1e-10){
                
                if (walker.location == Walker::intra){  

                    if (!isPosInsideAxon_(walker.pos_v, EPS_VAL)){
                        //cout << " is not inside axon " << id << " pos : " << walker.pos_v << endl;
                        collision.col_location = Collision::outside;
                    }
                    else{
                        //cout <<"tricky situation (intra)"<< endl;
                        continue;
                    }  
                }
                else{
                    collision.col_location = Collision::outside;
                } 
                
                collision.col_location = Collision::outside;
            }
            else {
                collision.col_location = Collision::unknown;
            }

            collision.type = Collision::hit;
            collision.collision_point = pos;
            collision.obstacle_ind = id;
            collision.obstacle_type = 0;
            collision.t = distance;
            collision.bounced_direction = -ray + 2.0 * normal * rn;
            collision.perm_crossing = 0.0;

            // Handle permeability

            if (this->percolation > 0.0) {
                static std::mt19937 gen_perm(std::random_device{}()); // Random number generator for permeability
                std::uniform_real_distribution<double> udist(0, 1);

                double dynamic_percolation = (collision.col_location == Collision::inside)
                                                 ? this->prob_cross_i_e
                                                 : this->prob_cross_e_i;

                double u = udist(gen_perm); 

                if (dynamic_percolation > u) {
                    collision.t += EPS_VAL;
                    collision.collision_point = walker.pos_v + collision.t * step;

                    count_perc_crossings++;
                    collision.perm_crossing = dynamic_percolation;
                    collision.bounced_direction = step;

                    return true;
                }
            }

            return true;
        }

    }

    collision.type = Collision::null;
    return false;
}


bool Axon::isWalkerInsideAxon(const Walker &walker, const double &distance_to_be_inside){
    
    Eigen::Vector3d O = walker.pos_v;
    bool isinside = isPosInsideAxon_(O, distance_to_be_inside);
    
    return isinside;
}
std::vector<int> Axon::checkAxisForCollision(const Eigen::Vector3d &position, const double &distance_to_be_inside, const int &axis){

    std::vector<int> spheres_id_to_check;

    for (int i = 0; i < spheres.size(); ++i) {
        const auto &sphere = spheres[i];

        // Check if the position is within the extended bounds of the sphere along the given axis
        double min_i = sphere.P[axis] - sphere.radius;
        double max_i = sphere.P[axis] + sphere.radius;

        if (position[axis] >= min_i - distance_to_be_inside && position[axis] <= max_i + distance_to_be_inside) {
            spheres_id_to_check.push_back(i);
        }
    }

    return spheres_id_to_check;
}

std::vector<int> Axon::findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3) {
    std::unordered_set<int> set1(vec1.begin(), vec1.end());
    std::unordered_set<int> set2(vec2.begin(), vec2.end());
    std::unordered_set<int> set3(vec3.begin(), vec3.end());

    std::vector<int> result;
    for (int val : set1) {
        if (set2.count(val) && set3.count(val)) {
            result.push_back(val);
        }
    }
    return result;
}
 bool Axon::isPosInsideAxon_(const Eigen::Vector3d &position, const double &distance_to_be_inside){
  
    if (!isNearAxon(position, distance_to_be_inside)) {
        return false; // Position is outside all bounding boxes
    }

    // Collect candidate spheres across all axes
    std::vector<std::vector<int>> spheres_id_to_check(3);
    for (int axis = 0; axis < 3; ++axis) {
        spheres_id_to_check[axis] = checkAxisForCollision(position, distance_to_be_inside, axis);
        if (spheres_id_to_check[axis].empty()) {
            return false; // No candidates for this axis
        }
    }

    // Find common spheres across all axes
    std::vector<int> spheres_to_check = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
    for (int sphere_id : spheres_to_check) {
        Sphere sphere_to_check = spheres[sphere_id];
        if (sphere_to_check.minDistance(position) < distance_to_be_inside) {
            return true;
        }
    }

    return false;
}

bool Axon::FindSphereinAxon(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int> &sph_ids){
  
    //cout << "isSphereInsideAxon_ : " << id << endl;
    sph_ids.clear();
    if(!isNearAxon(position, distance_to_be_inside)){ // if near axon
        return false;
    }

    //cout << "is near axon : " << id << endl;
    std::vector<std::vector<int>> spheres_id_to_check;
    for (auto axis = 0; axis < 3; ++axis) {
        spheres_id_to_check.push_back(checkAxisForCollision(position,distance_to_be_inside, axis)); // check for collision along 1 axis
        if (spheres_id_to_check[axis].size() == 0){
            return false;
        }
    }
    // find common ids in all 3 axes
    std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
    for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) {
        Sphere sphere_to_check = spheres[spheres_to_check_all_axes[i]];
        if (sphere_to_check.minDistance(position) < distance_to_be_inside){
            sph_ids.push_back(spheres_to_check_all_axes[i]);
            
        }
    }
    spheres_id_to_check.clear();
    spheres_to_check_all_axes.clear();
    
    if (sph_ids.size()>0){
        return true;
    }
    else{
        return false;
    }
    
}

bool Axon::isInsideOneBox(const Eigen::Vector3d& position, const double &distance_to_be_inside, const Box &box){

    // Expand box dimensions by distance_to_be_inside
    double x_min = box.x_min - distance_to_be_inside;
    double x_max = box.x_max + distance_to_be_inside;
    double y_min = box.y_min - distance_to_be_inside;
    double y_max = box.y_max + distance_to_be_inside;
    double z_min = box.z_min - distance_to_be_inside;
    double z_max = box.z_max + distance_to_be_inside;

    // Check if the position is inside the expanded box
    if (position[0] < x_min || position[0] > x_max) return false; // x-axis
    if (position[1] < y_min || position[1] > y_max) return false; // y-axis
    if (position[2] < z_min || position[2] > z_max) return false; // z-axis

    return true; // Inside all limits
}

bool Axon::isInsideBox(const Eigen::Vector3d& position, const double &distance_to_be_inside){
    for (auto box : boxes) {
        if (isInsideOneBox(position, distance_to_be_inside, box)) {
            return true;
        }
    }
    return false;
} 


bool Axon::isNearAxon(const Eigen::Vector3d &position, const double &distance_to_be_inside){

    if (isInsideBox(position, distance_to_be_inside)){
        return true;
    }
    
    return false;
}

bool Axon::isNearAxon(const Walker &walker, const double &distance_to_be_inside){

    Eigen::Vector3d position = walker.pos_v;

    return isNearAxon(position,distance_to_be_inside);
}

double Axon::minDistance(const Walker &w){
    //Origin of the ray
    Vector3d O = w.pos_v;
    return minDistance(O);
}

// Function to compute the distance from a point to a single box
double Axon::distanceToBox(const Eigen::Vector3d& O, const Box &box) {
    double dx = std::max({box.x_min - O[0], 0.0, O[0] - box.x_max});
    double dy = std::max({box.y_min - O[1], 0.0, O[1] - box.y_max});
    double dz = std::max({box.z_min - O[2], 0.0, O[2] - box.z_max});
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to compute the distance from a point to a single box
double Axon::distanceToBoxes(const Eigen::Vector3d& O) {
    double min_dist = INFINITY_VALUE;
    for (auto box : boxes) {
        double dist = distanceToBox(O, box);
        if (dist < min_dist) {
            min_dist = dist;
        }
    }
    return min_dist;
}


double Axon::minDistance(const Eigen::Vector3d &O) {

    double dist = distanceToBoxes(O);

    return dist;
}