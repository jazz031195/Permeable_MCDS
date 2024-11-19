#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include <numeric> // Add this line to include std::iota
#include <unordered_set>

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
}


void Axon::set_spheres(std::vector<Sphere> spheres_to_add){

    boxes.clear();

    int nbr_boxes = 1;
   
    int nbr_spheres_per_box = spheres_to_add.size() / nbr_boxes;

    // Check if there are spheres to add
    if (spheres_to_add.empty()) {
        std::cout << "No spheres to add." << std::endl;
        return;
    }

    // Create boxes
    for (int i = 0; i < nbr_boxes; i++) {
        Box box = Box();

        // Properly initialize box boundaries to 0.0
        box.x_min = 0.0;
        box.x_max = 0.0;
        box.y_min = 0.0;
        box.y_max = 0.0;
        box.z_min = 0.0;
        box.z_max = 0.0;

        for (int j = 0; j < nbr_spheres_per_box; ++j) {
            int index = j + nbr_spheres_per_box * i;  // Fix indexing calculation

            // Check for valid index
            if (index < spheres_to_add.size()) {
                spheres_to_add[index].object_id = id;
                spheres_to_add[index].id = index;

                if (j == 0) {
                    // Initialize box dimensions based on the first sphere
                    box.x_min = spheres_to_add[index].P[0] - spheres_to_add[index].radius;
                    box.x_max = spheres_to_add[index].P[0] + spheres_to_add[index].radius;
                    box.y_min = spheres_to_add[index].P[1] - spheres_to_add[index].radius;
                    box.y_max = spheres_to_add[index].P[1] + spheres_to_add[index].radius;
                    box.z_min = spheres_to_add[index].P[2] - spheres_to_add[index].radius;
                    box.z_max = spheres_to_add[index].P[2] + spheres_to_add[index].radius;
                } else {
                    // Adjust box dimensions for each subsequent sphere
                    if (spheres_to_add[index].P[0] + spheres_to_add[index].radius > box.x_max) {
                        box.x_max = spheres_to_add[index].P[0] + spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[0] - spheres_to_add[index].radius < box.x_min) {
                        box.x_min = spheres_to_add[index].P[0] - spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[1] + spheres_to_add[index].radius > box.y_max) {
                        box.y_max = spheres_to_add[index].P[1] + spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[1] - spheres_to_add[index].radius < box.y_min) {
                        box.y_min = spheres_to_add[index].P[1] - spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[2] + spheres_to_add[index].radius > box.z_max) {
                        box.z_max = spheres_to_add[index].P[2] + spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[2] - spheres_to_add[index].radius < box.z_min) {
                        box.z_min = spheres_to_add[index].P[2] - spheres_to_add[index].radius;
                    }
                }
            }
        }
        boxes.push_back(box);
    }

    // Handle leftover spheres
    int left_overs = spheres_to_add.size() % nbr_boxes;

    if (left_overs > 0) {
        
        Box last_box = Box();
        for (int j = 0; j < left_overs; ++j) {
            int index = j + nbr_boxes * nbr_spheres_per_box;  // Correct index calculation

            // Check for valid index
            if (index < spheres_to_add.size()) {
                spheres_to_add[index].object_id = id;
                spheres_to_add[index].id = index;

                if (j == 0) {
                    // Initialize last_box dimensions based on the first leftover sphere
                    last_box.x_min = spheres_to_add[index].P[0] - spheres_to_add[index].radius;
                    last_box.x_max = spheres_to_add[index].P[0] + spheres_to_add[index].radius;
                    last_box.y_min = spheres_to_add[index].P[1] - spheres_to_add[index].radius;
                    last_box.y_max = spheres_to_add[index].P[1] + spheres_to_add[index].radius;
                    last_box.z_min = spheres_to_add[index].P[2] - spheres_to_add[index].radius;
                    last_box.z_max = spheres_to_add[index].P[2] + spheres_to_add[index].radius;
                } else {
                    // Adjust last_box dimensions for each subsequent leftover sphere
                    if (spheres_to_add[index].P[0] + spheres_to_add[index].radius > last_box.x_max) {
                        last_box.x_max = spheres_to_add[index].P[0] + spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[0] - spheres_to_add[index].radius < last_box.x_min) {
                        last_box.x_min = spheres_to_add[index].P[0] - spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[1] + spheres_to_add[index].radius > last_box.y_max) {
                        last_box.y_max = spheres_to_add[index].P[1] + spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[1] - spheres_to_add[index].radius < last_box.y_min) {
                        last_box.y_min = spheres_to_add[index].P[1] - spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[2] + spheres_to_add[index].radius > last_box.z_max) {
                        last_box.z_max = spheres_to_add[index].P[2] + spheres_to_add[index].radius;
                    }
                    if (spheres_to_add[index].P[2] - spheres_to_add[index].radius < last_box.z_min) {
                        last_box.z_min = spheres_to_add[index].P[2] - spheres_to_add[index].radius;
                    }
                }
            }
        }
        boxes.push_back(last_box);
    }

    // Set the beginning and end positions
    this->begin = spheres_to_add[0].P;
    this->end = spheres_to_add[spheres_to_add.size() - 1].P;
    this->spheres = spheres_to_add;
    this->radius = spheres_to_add[0].radius;


}


bool Axon::intersection_sphere_vector(double &t1, double &t2, const Sphere &s, const Eigen::Vector3d &step, const Eigen::Vector3d &pos) {
    //https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    // Calculate vector from position to sphere center
    Eigen::Vector3d m = pos - s.P;
    double rad2 = s.radius * s.radius;

    // Quadratic equation coefficients
    double b = m.dot(step);
    double c = m.squaredNorm() - rad2;
    double discr = b * b - c; // a = 1 is implicit, omitted for optimization

    // If discriminant is negative, no intersection
    if (discr < 0.0) {
        return false;
    }

    double sqrt_discr = sqrt(discr);
    t1 = -b + sqrt_discr;
    t2 = -b - sqrt_discr;

    return true;
}


void Axon::set_prob_crossings(double step_length_pref){

    double prob_cross_i_e_, prob_cross_e_i_;
    double dse, dsi;

    if (this->percolation > 0.0){
        // for axon object
        dse = sqrt(step_length_pref*this->diffusivity_e);
        dsi = sqrt(step_length_pref*this->diffusivity_i);

        prob_cross_i_e_ = percolation * dsi * 2. / 3. / this->diffusivity_i;
        prob_cross_e_i_ = percolation * dse * 2. / 3. / this->diffusivity_e; 

        this->prob_cross_e_i = prob_cross_e_i_ / (1.+ 0.5 * (prob_cross_e_i_ + prob_cross_i_e_));
        this->prob_cross_i_e = prob_cross_i_e_ / (1.+ 0.5 * (prob_cross_e_i_ + prob_cross_i_e_));
    }
            
    // for all spheres
    for(unsigned i= 0 ; i < spheres.size();i++){
        spheres[i].percolation = this->percolation;

        if(spheres[i].percolation > 0.0){
            spheres[i].diffusivity_e = this->diffusivity_e;
            spheres[i].diffusivity_i = this->diffusivity_i;

            spheres[i].prob_cross_e_i = this->prob_cross_e_i;
            spheres[i].prob_cross_i_e = this->prob_cross_i_e;
            
        }
    }
}

void Axon::find_all_intersections(const Walker &walker, const Eigen::Vector3d &step, const double &step_length,
                                  std::vector<double> &dist_intersections, std::vector<int> &spheres_ids) {
    dist_intersections.clear();
    spheres_ids.clear();

    std::vector<int> sph_ids_walker_is_inside;
    Eigen::Vector3d pos = walker.pos_v;

    // Find spheres the walker is near
    if (!FindSphereinAxon(pos, 2*step_length, sph_ids_walker_is_inside)) {
        return;
    }

    // Process each sphere the walker is near
    for (int sphere_id : sph_ids_walker_is_inside) {
        Sphere sphere_to_check = spheres[sphere_id];

        // Skip if the sphere is too far
        if (sphere_to_check.minDistance(pos) > step_length) {
            continue;
        }

        double t1, t2;
        if (intersection_sphere_vector(t1, t2, sphere_to_check, step, pos)) {
            // Push valid intersections based on walker status
            if (walker.status == Walker::bouncing) {
                if (t1 > EPS_VAL) {
                    dist_intersections.push_back(t1);
                    spheres_ids.push_back(sphere_id);
                }
                if (t2 > EPS_VAL) {
                    dist_intersections.push_back(t2);
                    spheres_ids.push_back(sphere_id);
                }
            } else {
                if (t1 > 0) {
                    dist_intersections.push_back(t1);
                    spheres_ids.push_back(sphere_id);
                }
                if (t2 > 0) {
                    dist_intersections.push_back(t2);
                    spheres_ids.push_back(sphere_id);
                }
            }
        }
    }
}

bool Axon::checkCollision(Walker &walker, Eigen::Vector3d &step, const double &step_length, Collision &collision) {
    // Distances to intersections and corresponding sphere IDs
    std::vector<double> dist_intersections;
    std::vector<int> sphere_ids;

    walker.previous_location = walker.location;

    // Find all intersections
    find_all_intersections(walker, step, step_length + barrier_tickness, dist_intersections, sphere_ids);

    if (dist_intersections.empty()) {
        // Handle case with no intersections
        if (walker.location == Walker::intra && !isPosInsideAxon_(walker.pos_v, EPS_VAL)) {
            //cout << "Walker is outside axon" << endl;
            collision.col_location = Collision::outside;
            walker.in_obj_index = -1;
            walker.in_obj_type = -1;
            walker.location = Walker::extra;
            assert(0);
            return true;
        }

        collision.type = Collision::null;
        return false;
    }

    // Sort distances and corresponding sphere IDs
    std::vector<size_t> sorted_indices(dist_intersections.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(), [&dist_intersections](size_t i1, size_t i2) {
        return dist_intersections[i1] < dist_intersections[i2];
    });
    
    for (size_t i : sorted_indices) {
        
        double distance = dist_intersections[i];

        Eigen::Vector3d pos = walker.pos_v + distance * step;

        // Check if position is near the edge
        bool is_near_edge = (walker.location == Walker::intra) ? !isPosInsideAxon_(pos, -EPS_VAL) : true;

        if (is_near_edge && distance <= step_length + barrier_tickness) {

            Sphere sphere = spheres[sphere_ids[i]];

            collision.type = Collision::hit;
            collision.collision_point = pos;
            collision.obstacle_ind = id;
            collision.t = distance;

            // Compute normal and bounced direction
            Eigen::Vector3d normal = (collision.collision_point - sphere.P).normalized();
            Eigen::Vector3d ray = (-collision.t * step).normalized();
            double rn = ray.dot(normal);

            collision.bounced_direction = -ray + 2.0 * normal * rn;

            if (rn < -1e-10){
                
                if (walker.location == Walker::extra){ 
                    bool isinside_ = isPosInsideAxon_(walker.pos_v, -EPS_VAL);
                    if (isinside_){
                        collision.col_location = Collision::inside;
                        walker.in_obj_index = id;
                        walker.in_obj_type = 0;
                        walker.location = Walker::intra;
                    }
                    else{
                        collision.col_location = Collision::outside;
                        walker.in_obj_index = -1;
                        walker.in_obj_type = -1;
                        walker.location = Walker::extra;
                    }
                }
                else{
                    collision.col_location = Collision::inside;
                    walker.in_obj_index = id;
                    walker.in_obj_type = 0;
                    walker.location = Walker::intra;
                }
            
            }
            //outside
            else if (rn > 1e-10){
                if (walker.location == Walker::intra){ 
                    bool isinside_ = isPosInsideAxon_(walker.pos_v, EPS_VAL);
                    if (!isinside_){
                        collision.col_location = Collision::outside;
                        walker.in_obj_index = -1;
                        walker.in_obj_type = -1;
                        walker.location = Walker::extra;
                    }
                    else{
                        collision.col_location = Collision::inside;
                        walker.in_obj_index = id;
                        walker.in_obj_type = 0;
                        walker.location = Walker::intra;
                    }
                }
                else{
                    collision.col_location = Collision::outside;
                    walker.in_obj_index = -1;
                    walker.in_obj_type = -1;
                    walker.location = Walker::extra;
                }
            }
            else {
                collision.col_location = Collision::unknown;
            }

            // Handle permeability
            if (percolation > 0.0) {
                static std::mt19937 gen_perm(std::random_device{}());
                std::uniform_real_distribution<double> udist(0, 1);

                double dynamic_percolation = (collision.col_location == Collision::inside)
                                                 ? prob_cross_i_e
                                                 : prob_cross_e_i;

                if (dynamic_percolation > udist(gen_perm)) {
                    count_perc_crossings++;
                    collision.perm_crossing = dynamic_percolation;
                    collision.bounced_direction = step;
                    return true;
                }
            }

            collision.perm_crossing = 0.0;
            return true;
        }
    }

    collision.type = Collision::null;
    return false;
}


bool Axon::isWalkerInsideAxon(Walker &walker, double distance_to_be_inside){
    
    Eigen::Vector3d O;
    walker.getVoxelPosition(O);
    bool isinside = isPosInsideAxon_(O, distance_to_be_inside);
    
    return isinside;
}
std::vector<int> Axon::checkAxisForCollision(Eigen::Vector3d position, double distance_to_be_inside, int axis){

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
 bool Axon::isPosInsideAxon_(Eigen::Vector3d position, double distance_to_be_inside){
  
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
        if (sphere_to_check.minDistance(position) <= distance_to_be_inside) {
            return true;
        }
    }

    return false;
}

bool Axon::FindSphereinAxon(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int> &sph_ids){
  
    //cout << "isSphereInsideAxon_ : " << id << endl;
    sph_ids.clear();
    if(isNearAxon(position, distance_to_be_inside)){ // if near axon
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
            if (sphere_to_check.minDistance(position) <= distance_to_be_inside){
                sph_ids.push_back(spheres_to_check_all_axes[i]);
                
            }
        }
        spheres_id_to_check.clear();
        spheres_to_check_all_axes.clear();
    }
    if (sph_ids.size()>0){
        return true;
    }
    else{
        return false;
    }
    
}

bool Axon::isInsideBox(const int& i, const Eigen::Vector3d& position, const double &distance_to_be_inside){

    // Expand box dimensions by distance_to_be_inside
    double x_min = boxes[i].x_min - distance_to_be_inside;
    double x_max = boxes[i].x_max + distance_to_be_inside;
    double y_min = boxes[i].y_min - distance_to_be_inside;
    double y_max = boxes[i].y_max + distance_to_be_inside;
    double z_min = boxes[i].z_min - distance_to_be_inside;
    double z_max = boxes[i].z_max + distance_to_be_inside;

    // Check if the position is inside the expanded box
    if (position[0] < x_min || position[0] > x_max) return false; // x-axis
    if (position[1] < y_min || position[1] > y_max) return false; // y-axis
    if (position[2] < z_min || position[2] > z_max) return false; // z-axis

    return true; // Inside all limits
}

bool Axon::isNearAxon(Eigen::Vector3d position, double distance_to_be_inside){
    
    for (int i = 0; i < boxes.size(); ++i) {
        if (isInsideBox(i, position, distance_to_be_inside)){
            return true;
        }
    }

    return false;
}

bool Axon::isNearAxon(Walker walker, double distance_to_be_inside){

    Eigen::Vector3d position;
    walker.getVoxelPosition(position);

    return isNearAxon(position,distance_to_be_inside);

}


double Axon::minDistance(Walker &w){
    //Origin of the ray
    Vector3d O;
    w.getVoxelPosition(O);

    return minDistance(O);
}

// Function to compute the distance from a point to a single box
double Axon::distanceToBox(const int& i, const Eigen::Vector3d& O) {
    double dx = std::max({boxes[i].x_min - O[0], 0.0, O[0] - boxes[i].x_max});
    double dy = std::max({boxes[i].y_min - O[1], 0.0, O[1] - boxes[i].y_max});
    double dz = std::max({boxes[i].z_min - O[2], 0.0, O[2] - boxes[i].z_max});
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Axon::minDistance(const Eigen::Vector3d &O) {
    // Check if there are no boxes
    if (boxes.empty()) {
        std::cerr << "Error: No boxes available to calculate distance." << std::endl;
        return std::numeric_limits<double>::infinity();  // Return a large value to indicate no boxes
    }

    double min_distance = std::numeric_limits<double>::max();
    
    for (int i = 0; i < boxes.size(); ++i) {
        double dist = distanceToBox(i, O);
        if (dist < min_distance) {
            min_distance = dist;
        }
    }
    
    return min_distance;
}