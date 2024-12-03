#include "Glial.h"
#include "Eigen/Dense"
#include "constants.h"
#include <numeric>
#include <unordered_map>
#include <map>
#include <unordered_set>

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
    boxes = gl.boxes;
    big_box = gl.big_box;

    percolation = gl.percolation;
    prob_cross_e_i = gl.prob_cross_e_i;
    prob_cross_i_e = gl.prob_cross_i_e;
    diffusivity_i = gl.diffusivity_i;
    diffusivity_e = gl.diffusivity_e;
    count_perc_crossings = gl.count_perc_crossings;

};

void Glial::set_spheres(std::vector<Sphere> &spheres_to_add) {
    // Clear existing boxes and initialize variables
    boxes.clear();
    processes.clear();

    big_box = {
        soma.P[0] - soma.radius, soma.P[0] + soma.radius,
        soma.P[1] - soma.radius, soma.P[1] + soma.radius,
        soma.P[2] - soma.radius, soma.P[2] + soma.radius
    };

    boxes.push_back(big_box);

    std::map<int, Box> branch_boxes; // Map branch_id to Box

    for (const auto &sphere : spheres_to_add) {

        double x_min = sphere.P[0] - sphere.radius;
        double x_max = sphere.P[0] + sphere.radius;
        double y_min = sphere.P[1] - sphere.radius;
        double y_max = sphere.P[1] + sphere.radius;
        double z_min = sphere.P[2] - sphere.radius;
        double z_max = sphere.P[2] + sphere.radius;


        if (x_min< big_box.x_min) big_box.x_min = x_min;
        if (x_max> big_box.x_max) big_box.x_max = x_max;
        if (y_min< big_box.y_min) big_box.y_min = y_min;
        if (y_max> big_box.y_max) big_box.y_max = y_max;
        if (z_min< big_box.z_min) big_box.z_min = z_min;
        if (z_max> big_box.z_max) big_box.z_max = z_max;

        // Handle processes
        processes.push_back(sphere);

        // Update or create the bounding box for this branch_id
        if (branch_boxes.find(sphere.branch_id) == branch_boxes.end()) {

            // Create new box
            branch_boxes[sphere.branch_id] = {
                x_min, x_max, y_min, y_max, z_min, z_max
            };

        } else {
            // Update existing box
            Box &box = branch_boxes[sphere.branch_id];
            if (x_min < box.x_min) box.x_min = x_min;
            if (x_max > box.x_max) box.x_max = x_max;
            if (y_min < box.y_min) box.y_min = y_min;
            if (y_max > box.y_max) box.y_max = y_max;
            if (z_min < box.z_min) box.z_min = z_min;
            if (z_max > box.z_max) box.z_max = z_max;

        }
        
    }

    // Add all branch boxes to the final list of boxes
    for (const auto &entry : branch_boxes) {
        boxes.push_back(entry.second);
    }

}

bool Glial::isInsideBox(const int &i, const Eigen::Vector3d &position, const double &distance_to_be_inside) {
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

bool Glial::isInsideBigBox(const Eigen::Vector3d &position, const double &distance_to_be_inside) {
    // Expand box dimensions by distance_to_be_inside
    double x_min = big_box.x_min - distance_to_be_inside;
    double x_max = big_box.x_max + distance_to_be_inside;
    double y_min = big_box.y_min - distance_to_be_inside;
    double y_max = big_box.y_max + distance_to_be_inside;
    double z_min = big_box.z_min - distance_to_be_inside;
    double z_max = big_box.z_max + distance_to_be_inside;

    // Check if the position is inside the expanded box
    if (position[0] < x_min || position[0] > x_max) return false; // x-axis
    if (position[1] < y_min || position[1] > y_max) return false; // y-axis
    if (position[2] < z_min || position[2] > z_max) return false; // z-axis

    return true; // Inside all limits
}


bool Glial::isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int> &branches) {

    branches.clear();
    if (isInsideBigBox(position, distance_to_be_inside)) {
        for (int i = 0; i < boxes.size(); ++i) {
            if (isInsideBox(i, position, distance_to_be_inside)) {
                branches.push_back(i-1);
                //return true; // Near at least one box
            }
        }
    }

    if (!branches.empty()) return true;

    return false; // Not near any box
}

std::vector<int> Glial::checkAxisForCollision(const Eigen::Vector3d &position, double distance_to_be_inside, int axis, const std::vector<int> &branches) {
    std::vector<int> spheres_id_to_check;

    // Collect relevant spheres
    std::vector<Sphere> relevant_spheres = {soma};
    relevant_spheres.insert(relevant_spheres.end(), processes.begin(), processes.end());

    for (int i = 0; i < relevant_spheres.size(); ++i) {
        const auto &sphere = relevant_spheres[i];

        // Check if the sphere's branch_id is in the provided branches
        if (std::find(branches.begin(), branches.end(), sphere.branch_id) != branches.end()) {
            
            // Check if the position is within the extended bounds of the sphere along the given axis
            double min_i = sphere.P[axis] - sphere.radius;
            double max_i = sphere.P[axis] + sphere.radius;

            if (position[axis] >= min_i - distance_to_be_inside && position[axis] <= max_i + distance_to_be_inside) {
                spheres_id_to_check.push_back(i);
            }
        }
        
    }

    return spheres_id_to_check;
}


std::vector<int> Glial::findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3) {
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

bool Glial::isPosInsideGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside) {
    
    std::vector<int> branches;
    if (!isNearGlialCell(position, distance_to_be_inside, branches)) {
        return false; // Position is outside all bounding boxes
    }

    // Collect candidate spheres across all axes
    std::vector<std::vector<int>> spheres_id_to_check(3);
    for (int axis = 0; axis < 3; ++axis) {
        spheres_id_to_check[axis] = checkAxisForCollision(position, distance_to_be_inside, axis, branches);
        if (spheres_id_to_check[axis].empty()) {
            return false; // No candidates for this axis
        }
    }

    // Find common spheres across all axes
    std::vector<int> spheres_to_check = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
    for (int sphere_id : spheres_to_check) {
        Sphere sphere_to_check = (sphere_id == 0) ? soma : processes[sphere_id - 1];
        if (sphere_to_check.minDistance(position) <= distance_to_be_inside) {
            return true;
        }
    }

    return false;
}

bool Glial::FindSphereinGlial(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int> &sph_ids) {
    sph_ids.clear();

    // Determine nearby branches
    std::vector<int> branches;
    if (!isNearGlialCell(position, distance_to_be_inside, branches)) {
        return false;
    }
    

    // Collect potential sphere IDs from all axes
    std::vector<std::vector<int>> spheres_id_to_check(3);
    for (int axis = 0; axis < 3; ++axis) {
        spheres_id_to_check[axis] = checkAxisForCollision(position, distance_to_be_inside, axis, branches);
        if (spheres_id_to_check[axis].empty()) {
            return false; // No collisions along this axis
        }
    }

    // Find common sphere IDs across all axes
    std::vector<int> spheres_to_check_all_axes = findCommonIntegers(
        spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]
    );

    // Check each common sphere for proximity
    for (int sphere_id : spheres_to_check_all_axes) {
        Sphere sphere_to_check = (sphere_id == 0) ? soma : processes[sphere_id - 1];
        if (sphere_to_check.minDistance(position) <= distance_to_be_inside) {
            sph_ids.push_back(sphere_id);

        }
    }

    return !sph_ids.empty();
}


bool Glial::intersection_sphere_vector(double &t1, double &t2, const Sphere &s, const Eigen::Vector3d &step, const Eigen::Vector3d &pos) {
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

void Glial::find_all_intersections(const Walker &walker, const Eigen::Vector3d &step, const double &distance,
                                   std::vector<double> &dist_intersections, std::vector<int> &spheres_ids) {
    dist_intersections.clear();
    spheres_ids.clear();

    std::vector<int> sph_ids_walker_is_inside;
    Eigen::Vector3d pos = walker.pos_v;

    // Find spheres the walker is near
    if (!FindSphereinGlial(pos, distance, sph_ids_walker_is_inside)) {
        return;
    }

    // Process each sphere the walker is near
    for (int sphere_id : sph_ids_walker_is_inside) {
        Sphere sphere_to_check = (sphere_id == 0) ? soma : processes[sphere_id - 1];

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

bool Glial::checkCollision(const Walker &walker, Eigen::Vector3d &step, const double &step_length, Collision &collision) {
    // Distances to intersections and corresponding sphere IDs
    std::vector<double> dist_intersections;
    std::vector<int> sphere_ids;


    // Find all intersections
    find_all_intersections(walker, step, step_length + barrier_tickness, dist_intersections, sphere_ids);

    if (dist_intersections.empty()) {
        // Handle case with no intersections
        if (walker.location == Walker::intra && !isPosInsideGlialCell(walker.pos_v, EPS_VAL)) {
            //cout << "Walker is outside axon" << endl;
            collision.type = Collision::hit;
            collision.col_location = Collision::outside;
            collision.collision_point = walker.pos_v;
            collision.t = 1e-9;
            collision.perm_crossing = 0.0;
            collision.bounced_direction = step;
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
        //out << "---" << endl;
        double distance = dist_intersections[i];
        //cout << "distance : " << distance << endl;
        Eigen::Vector3d pos = walker.pos_v + distance * step;

        // Check if position is near the edge
        bool is_near_edge = (walker.location == Walker::intra) ? !isPosInsideGlialCell(pos, -EPS_VAL) : true;

        if (is_near_edge && distance <= step_length + barrier_tickness) {
            //cout << "is near edge" << endl;
            Sphere sphere = (sphere_ids[i] == 0) ? soma : processes[sphere_ids[i] - 1];

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
                    if (!isPosInsideGlialCell(walker.pos_v, EPS_VAL)){
                        collision.col_location = Collision::outside;
                    }
                    else{
                        continue;
                    }  
                }
                else{
                    collision.col_location = Collision::outside;
                } 
            }
            else {
                collision.col_location = Collision::unknown;
            }

            collision.type = Collision::hit;
            collision.collision_point = pos;
            collision.obstacle_ind = id;
            collision.obstacle_type = 1;
            collision.t = distance;
            collision.bounced_direction = -ray + 2.0 * normal * rn;


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
        processes[i].percolation = this->percolation;

        if(processes[i].percolation > 0.0){
            processes[i].prob_cross_e_i = this->prob_cross_e_i;
            processes[i].prob_cross_i_e = this->prob_cross_i_e;

            processes[i].diffusivity_e = this->diffusivity_e;
            processes[i].diffusivity_i = this->diffusivity_i;
            
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

// Function to compute the distance from a point to a single box
double Glial::distanceToBox(const int& i, const Eigen::Vector3d& O) {
    double dx = std::max({boxes[i].x_min - O[0], 0.0, O[0] - boxes[i].x_max});
    double dy = std::max({boxes[i].y_min - O[1], 0.0, O[1] - boxes[i].y_max});
    double dz = std::max({boxes[i].z_min - O[2], 0.0, O[2] - boxes[i].z_max});
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}


double Glial::distanceToBigBox(const Eigen::Vector3d& O) {
    double dx = std::max({big_box.x_min - O[0], 0.0, O[0] - big_box.x_max});
    double dy = std::max({big_box.y_min - O[1], 0.0, O[1] - big_box.y_max});
    double dz = std::max({big_box.z_min - O[2], 0.0, O[2] - big_box.z_max});
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

double Glial::minDistance(const Walker &w){

    // Check if there are no boxes
    if (boxes.empty()) {
        std::cerr << "Error: No boxes available to calculate distance." << std::endl;
        return std::numeric_limits<double>::infinity();  // Return a large value to indicate no boxes
    }

    // Get the walker's position
    Vector3d O = w.pos_v;

    double min_distance = std::numeric_limits<double>::max();
    
    // Check distance to big box
    //double min_distance = distanceToBigBox(O);

    for (int i = 0; i < boxes.size(); ++i) {
        double dist = distanceToBox(i, O);
        if (dist < min_distance) {
            min_distance = dist;
        }
    }
    
    return min_distance;
}


