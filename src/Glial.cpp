#include "Glial.h"
#include "Eigen/Dense"
#include "constants.h"
#include <numeric>
#include <unordered_map>

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

};

void Glial::set_spheres(const std::vector<Sphere> &spheres_to_add) {
    // Clear existing boxes and initialize variables
    boxes.clear();
    soma = Sphere();
    processes.clear();

    std::unordered_map<int, Box> branch_boxes; // Map branch_id to Box

    for (const auto &sphere : spheres_to_add) {
        if (sphere.branch_id == -1) {
            // Handle soma
            soma = sphere;
            Box soma_box = {
                sphere.P[0] - sphere.radius, sphere.P[0] + sphere.radius,
                sphere.P[1] - sphere.radius, sphere.P[1] + sphere.radius,
                sphere.P[2] - sphere.radius, sphere.P[2] + sphere.radius
            };
            boxes.push_back(soma_box);
        } else {
            // Handle processes
            processes.push_back(sphere);

            // Update or create the bounding box for this branch_id
            if (branch_boxes.find(sphere.branch_id) == branch_boxes.end()) {
                // Create new box
                branch_boxes[sphere.branch_id] = {
                    sphere.P[0] - sphere.radius, sphere.P[0] + sphere.radius,
                    sphere.P[1] - sphere.radius, sphere.P[1] + sphere.radius,
                    sphere.P[2] - sphere.radius, sphere.P[2] + sphere.radius
                };
            } else {
                // Update existing box
                Box &box = branch_boxes[sphere.branch_id];
                box.x_min = std::min(box.x_min, sphere.P[0] - sphere.radius);
                box.x_max = std::max(box.x_max, sphere.P[0] + sphere.radius);
                box.y_min = std::min(box.y_min, sphere.P[1] - sphere.radius);
                box.y_max = std::max(box.y_max, sphere.P[1] + sphere.radius);
                box.z_min = std::min(box.z_min, sphere.P[2] - sphere.radius);
                box.z_max = std::max(box.z_max, sphere.P[2] + sphere.radius);
            }
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

bool Glial::isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside) {
    // Check if the position is near any box
    for (int i = 0; i < boxes.size(); ++i) {
        if (isInsideBox(i, position, distance_to_be_inside)) {
            return true;
        }
    }

    return false; // Not near any box
}

std::vector<int> Glial::checkAxisForCollision(Eigen::Vector3d position, double distance_to_be_inside, int axis){

	std::vector<int> spheres_id_to_check;
    std::vector<Sphere> all_spheres = {soma};
    all_spheres.insert(all_spheres.end(), processes.begin(), processes.end());

    //cout << "all_spheres.size() : " << all_spheres.size() << endl;
	for (auto i = 0; i < all_spheres.size(); ++i) {

            double min_i = all_spheres[i].P[axis] - all_spheres[i].radius;

			if (min_i - distance_to_be_inside > position[axis]) {
				continue;
			}
			else {
				double max_i = all_spheres[i].P[axis] + all_spheres[i].radius;

                if (position[axis]> max_i + distance_to_be_inside) {
                    continue;
                }
                else{
                    spheres_id_to_check.push_back(i);
                }
			}
	}

    return spheres_id_to_check;
}
std::vector<int> Glial::findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3) {
    std::vector<int> result;
    std::set_intersection(vec1.begin(), vec1.end(),
                          vec2.begin(), vec2.end(),
                          std::back_inserter(result));
    std::vector<int> commonIntegers;
    std::set_intersection(result.begin(), result.end(),
                          vec3.begin(), vec3.end(),
                          std::back_inserter(commonIntegers));
    return commonIntegers;
}

bool Glial::isPosInsideGlialCell(const Eigen::Vector3d& position, const double& distance_to_be_inside){

    if (isNearGlialCell(position, distance_to_be_inside)){
  
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
            if (spheres_to_check_all_axes[i] == 0){
                // soma
                Sphere sphere_to_check = soma;
                if (sphere_to_check.minDistance(position) <= distance_to_be_inside){
                    return true;
                }

            }
            else{
                //processes
                Sphere sphere_to_check = processes[spheres_to_check_all_axes[i]-1];
                if (sphere_to_check.minDistance(position) <= distance_to_be_inside){
                    return true;
                }
            }
        }
        spheres_id_to_check.clear();
        spheres_to_check_all_axes.clear();
    }

    return false;
}


bool Glial::FindSphereinGlial(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int> &sph_ids){
  
    //cout << "isSphereInsideAxon_ : " << id << endl;
    sph_ids.clear();

    if (isNearGlialCell(position, distance_to_be_inside)){
  
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
            if (spheres_to_check_all_axes[i] == 0){
                // soma
                Sphere sphere_to_check = soma;
                if (sphere_to_check.minDistance(position) <= distance_to_be_inside){
                    sph_ids.push_back(spheres_to_check_all_axes[i] );
                }

            }
            else{
                //processes
                Sphere sphere_to_check = processes[spheres_to_check_all_axes[i]-1];
                if (sphere_to_check.minDistance(position) <= distance_to_be_inside){
                    sph_ids.push_back(spheres_to_check_all_axes[i]);
                }
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



bool Glial::intersection_sphere_vector(double &t1, double &t2, Sphere &s, Eigen::Vector3d &step, const double &step_length, const Eigen::Vector3d &pos){
    //https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    Eigen::Vector3d m = pos - s.P;
    double rad = s.radius;

    double a = 1;
    double b = (m.dot(step));
    double c = m.dot(m) - rad*rad;
    double discr = b*b - a*c;

    if (discr < 0.0 ){
        return false;
    }

    t1 = (-b + sqrt(discr))/(a);
    t2 = (-b - sqrt(discr))/(a);

    return true;

}

void Glial::find_all_intersections(const Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, std::vector<double>& dist_intersections, std::vector<int>& spheres_ids){

    dist_intersections.clear();
    spheres_ids.clear();

    std::vector<int> sph_ids_walker_is_inside;
    Eigen::Vector3d pos = walker.pos_v;
    bool isinside = FindSphereinGlial(pos, 2*step_lenght,sph_ids_walker_is_inside);

    if (!isinside){
        return;
    }

    for (auto i = 0; i < sph_ids_walker_is_inside.size() ; ++i) {


        Sphere sphere_to_check;

        if (sph_ids_walker_is_inside[i] == 0){
            //soma
            sphere_to_check = soma;
        }
        else{
            //processes
            sphere_to_check = processes[sph_ids_walker_is_inside[i]-1];

        }

        if (sphere_to_check.minDistance(walker.pos_v) > (step_lenght)){
            continue;
        }

        double t1, t2;
        bool intersect = intersection_sphere_vector(t1, t2, sphere_to_check, step, step_lenght, walker.pos_v);
        if (intersect){
            if (walker.status == Walker::bouncing){
                if (t1 > EPS_VAL){
                    dist_intersections.push_back(t1);
                    spheres_ids.push_back(sph_ids_walker_is_inside[i]);

                }
                if (t2 > EPS_VAL){
                    dist_intersections.push_back(t2);
                    spheres_ids.push_back(sph_ids_walker_is_inside[i]);
                }
            }
            else{
                if (t1 > 0){
                    dist_intersections.push_back(t1);
                    spheres_ids.push_back(sph_ids_walker_is_inside[i]);
                }
                if (t2 > 0){
                    dist_intersections.push_back(t2);
                    spheres_ids.push_back(sph_ids_walker_is_inside[i]);
                }
            }
        }
    }

    
}
bool Glial::checkCollision(Walker &walker, Eigen::Vector3d &step, const double &step_length, Collision &collision) {
    // Distances to intersections and corresponding sphere IDs
    std::vector<double> dist_intersections;
    std::vector<int> sphere_ids;

    walker.previous_location = walker.location;
    //cout << "**********" << endl;
    //cout << "is inside : " << isPosInsideGlialCell(walker.pos_v, EPS_VAL) << endl;
    
    // Find all intersections
    find_all_intersections(walker, step, step_length + barrier_tickness, dist_intersections, sphere_ids);

    if (dist_intersections.empty()) {
        // Handle case with no intersections
        if (walker.location == Walker::intra) {
            if (!isPosInsideGlialCell(walker.pos_v, EPS_VAL)) {
                collision.col_location = Collision::outside;
                walker.in_obj_index = -1;
                walker.in_obj_type = -1;
                walker.location = Walker::extra;
                cout << "outside" << endl;
                assert(false);
                return true;
            }
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

            collision.type = Collision::hit;
            collision.colision_point = pos;
            collision.obstacle_ind = id;
            collision.t = distance;

            // Compute normal and bounced direction
            Eigen::Vector3d normal = (collision.colision_point - sphere.P).normalized();
            Eigen::Vector3d ray = (-collision.t * step).normalized();
            double rn = ray.dot(normal);

            collision.bounced_direction = -ray + 2.0 * normal * rn;

            if (rn < -1e-10){
                if (walker.location == Walker::extra){ 
                    bool isinside_ = isPosInsideGlialCell(walker.pos_v, -EPS_VAL);
                    if (isinside_){
                        colision.col_location = Collision::inside;
                        walker.in_obj_index = id;
                        walker.in_obj_type = 0;
                        walker.location = Walker::intra;

                    }
                    else{
                        colision.col_location = Collision::outside;
                        walker.in_obj_index = -1;
                        walker.in_obj_type = -1;
                        walker.location = Walker::extra;
                    }
                }
                else{
                    colision.col_location = Collision::inside;
                    walker.in_obj_index = id;
                    walker.in_obj_type = 0;
                    walker.location = Walker::intra;
                }
                //cout << "rn :" << rn << endl;
                //assert(0);
                
            }
            //outside
            else if (rn > 1e-10){
                
                if (walker.location == Walker::intra){ 
                    bool isinside_ = isPosInsideGlialCell(walker.pos_v, EPS_VAL);
                    if (!isinside_){
                        colision.col_location = Collision::outside;
                        walker.in_obj_index = -1;
                        walker.in_obj_type = -1;
                        walker.location = Walker::extra;

                    }
                    else{
                        colision.col_location = Collision::inside;
                        walker.in_obj_index = id;
                        walker.in_obj_type = 0;
                        walker.location = Walker::intra;
                    }
                }
                else{
                    colision.col_location = Collision::outside;
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
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Glial::minDistance(Walker &w){

    // Check if there are no boxes
    if (boxes.empty()) {
        std::cerr << "Error: No boxes available to calculate distance." << std::endl;
        return std::numeric_limits<double>::infinity();  // Return a large value to indicate no boxes
    }

    // Get the walker's position
    Vector3d O = w.pos_v;

    double min_distance = std::numeric_limits<double>::max();
    
    for (int i = 0; i < boxes.size(); ++i) {
        double dist = distanceToBox(i, O);
        if (dist < min_distance) {
            min_distance = dist;
        }
    }
    
    return min_distance;
}
