#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>
#include <numeric> // Add this line to include std::iota


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

    int nbr_boxes = 10;

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


bool Axon::intersection_sphere_vector(double &t1, double &t2, Sphere &s, Eigen::Vector3d &step, const double &step_length, Eigen::Vector3d &pos){
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

double find_rn(double t, Eigen::Vector3d O, Eigen::Vector3d step, Eigen::Vector3d sphere_center){
    // center of sphere to collision point
    Eigen::Vector3d normal = ((t*step+O)- sphere_center).normalized();
    // - vector of step
    Eigen::Vector3d ray =  (-t*step).normalized();
    // projection onto normal
    double rn = ray.dot(normal);

    return rn;
}


void Axon::bouncing(Walker &walker, Collision &colision, double dist_to_collision, double step_lenght, int sphere_ind, Eigen::Vector3d step, int index_, std::vector<double> rns){
    if (rns[index_] <0){
        colision.col_location = Collision::inside;
        }
    else {
        colision.col_location = Collision::outside;          
    }

    colision.type = Collision::hit;
    colision.rn = rns[index_];  
    colision.obstacle_ind = id;
    colision.t = fmin(dist_to_collision,step_lenght);
    colision.colision_point = walker.pos_v + colision.t*step;
    walker.is_allowed_to_cross = false;

    // Membrane permeability    
    if((spheres[sphere_ind].percolation>0.0)){


        std::mt19937 gen_perm;          // Random engine for permeability
        std::random_device rd;
        gen_perm.seed(rd());
        std::uniform_real_distribution<double> udist(0,1);
                    
        double _percolation_ = udist(gen_perm); 

        double dynamic_percolation = 0.0;
                    
        if (colision.col_location == Collision::inside){ 
            dynamic_percolation =  spheres[sphere_ind].prob_cross_i_e; 
        } 

        else if (colision.col_location == Collision::outside){
            dynamic_percolation = spheres[sphere_ind].prob_cross_e_i;
        } 

        if( dynamic_percolation - _percolation_ > EPS_VAL ){            
            count_perc_crossings++;
            colision.perm_crossing      = _percolation_;
            colision.bounced_direction  = step; 
            walker.is_allowed_to_cross = true;

            // if crosses from extra to intra, save id of axon
            if (walker.location == Walker::extra){
                walker.in_obj_index = this->id;
                walker.in_obj_type = 0;
                        
            }
            
            return;
        }
                
    }
            
    colision.perm_crossing = 0.;

    /* For a sphere, normal direction is equal to colision point */
    //Normal point
    Eigen::Vector3d normal = (colision.colision_point- spheres[sphere_ind].P).normalized();
    Eigen::Vector3d temp_step = step;
    elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);
    colision.bounced_direction = temp_step.normalized();


} 

void Axon::find_all_intersections(const Walker &walker,  Eigen::Vector3d &step, const double &step_lenght, std::vector<double>& dist_intersections, std::vector<int>& spheres_ids){

    dist_intersections.clear();
    spheres_ids.clear();

    //new_sph_id_to_check.push_back(0);
    std::vector <int> spheres_walker_is_inside;
    bool is_intra = FindSphereinAxon(walker.pos_v, 2*step_lenght, spheres_walker_is_inside);

    if (!is_intra){
        return;
    }

    Eigen::Vector3d pos = walker.pos_v;
    if (spheres.size() == 0){
        cout << "spheres.size() == 0" << endl;
        assert(0);
    }

    for (auto i = 0; i < spheres_walker_is_inside.size(); ++i) {
    
        Sphere sphere_to_check = spheres[spheres_walker_is_inside[i]];

        double t1, t2;
        bool intersect = intersection_sphere_vector(t1, t2, sphere_to_check, step, step_lenght, pos);
        if (intersect){
            if (walker.status == Walker::bouncing){
                if (t1 > EPS_VAL){
                dist_intersections.push_back(t1);
                spheres_ids.push_back(spheres_walker_is_inside[i]);
                }
                if (t2 > EPS_VAL){
                    dist_intersections.push_back(t2);
                    spheres_ids.push_back(spheres_walker_is_inside[i]);
                }
            }
            else{
                if (t1 > 0){
                    dist_intersections.push_back(t1);
                    spheres_ids.push_back(spheres_walker_is_inside[i]);
                }
                if (t2 > 0){
                    dist_intersections.push_back(t2);
                    spheres_ids.push_back(spheres_walker_is_inside[i]);
                }
            }
        }
    }
}
bool Axon::checkCollision(Walker &walker,  Eigen::Vector3d &step, double step_lenght, Collision &colision){
    // distances to intersections
    std::vector<double> dist_intersections;

    walker.previous_location = walker.location;

    std::vector<int> spheres_ids;
    
    find_all_intersections(walker, step, step_lenght+barrier_tickness, dist_intersections, spheres_ids);

    if (dist_intersections.empty()){
        
        if (walker.location == Walker::intra){ 
            bool isinside = isPosInsideAxon_(walker.pos_v, EPS_VAL);
            if (!isinside){
                colision.col_location = Collision::outside;
                walker.in_obj_index = -1;
                walker.in_obj_type = -1;
                walker.location = Walker::extra;
                return true;
            }
        }
        
        colision.type = Collision::null;
        return false;
    }
    else{

        std::vector<int> idx(dist_intersections.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&dist_intersections](int i1, int i2) {return dist_intersections[i1] < dist_intersections[i2];});
        


        for (auto i = 0; i < idx.size(); ++i) {
            int index = idx[i];
            Eigen::Vector3d pos = walker.pos_v + dist_intersections[index]*step;

            bool isnearEdge = true;
            
            if (walker.location == Walker::intra){ 
                isnearEdge = !isPosInsideAxon_(pos, -barrier_tickness);
            }

            if (isnearEdge ){
                
                if (dist_intersections[index] <= step_lenght+barrier_tickness){
                    //cout <<"Chosen: step length :" <<step_lenght <<" i : " <<i <<" index : " << index << " dist_intersections[i] : " << dist_intersections[index] << " id_to_check[i] : " <<  spheres_ids[index] << " isnearEdge : " << isnearEdge << endl;

                    Sphere sph = spheres[spheres_ids[index]];
                    
                    //cout << "distance to sphere wall :" << sph.minDistance(pos)  << endl;
                    colision.type = Collision::hit;
                    colision.colision_point = pos;
                    colision.obstacle_ind = id;
                    colision.t = dist_intersections[index];


                    //Normal point
                    Eigen::Vector3d normal = (colision.colision_point - sph.P).normalized();
                    // Elastic bounce
                    Eigen::Vector3d temp_step = step;
                    Eigen::Vector3d ray =  (-colision.t*temp_step).normalized();//
                    double rn = ray.dot(normal);
                    temp_step = -ray + 2.0*normal*rn;
                    colision.perm_crossing = 0.;
                    //cout << "bounced_step : " << temp_step << " step : "<< step << endl;
                    //cout << "next_pos_ : " << walker.pos_v + colision.t*step << endl;
                    //cout << " collision point : " << colision.colision_point <<endl;

                    colision.bounced_direction = temp_step.normalized();
                    colision.rn = rn;

                    //cout << "rn :" << rn << endl;

                    //inside
                    if (rn < -1e-10){
                        if (walker.location == Walker::extra){ 
                            bool isinside_ = isPosInsideAxon_(walker.pos_v, -EPS_VAL);
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
                            bool isinside_ = isPosInsideAxon_(walker.pos_v, EPS_VAL);
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
                    else{
                        colision.col_location = Collision::unknown;
                    }

                    // Membrane permeability    
                    if((this->percolation>0.0)){
                        if(colision.type == Collision::hit && colision.col_location != Collision::voxel){

                            std::mt19937 gen_perm;          // Random engine for permeability
                            std::random_device rd;
                            gen_perm.seed(rd());
                            std::uniform_real_distribution<double> udist(0,1);
                            // generate a random number between 0 and 1
                            double _percolation_ = udist(gen_perm); 

                            double dynamic_percolation = 0.0;
                            
                            if (colision.col_location == Collision::inside){ 
                                dynamic_percolation =  this->prob_cross_i_e; 
                            } 

                            else if (colision.col_location == Collision::outside){
                                dynamic_percolation = this->prob_cross_e_i;
                            } 

                            if( dynamic_percolation - _percolation_ > EPS_VAL ){            
                                count_perc_crossings++;
                                colision.perm_crossing      = _percolation_;
                                colision.bounced_direction  = step; 
                                return true;
                            }
                        }  
                    }


                    colision.perm_crossing = 0.;
                    //cout <<"hit" << endl;

                    return true;
                }
                else{
                    //cout << "dist_intersections[index] > step_lenght+barrier_tickness" << endl;
                    colision.type = Collision::null;
                    return false;
                }
                
            }
            
        }
        colision.type = Collision::null;
        return false;
        
    }
    

}

bool Axon::isWalkerInsideAxon(Walker &walker, double distance_to_be_inside){
    
    Eigen::Vector3d O;
    walker.getVoxelPosition(O);
    bool isinside = isPosInsideAxon_(O, distance_to_be_inside);
    
    return isinside;
}


std::vector<int> Axon::checkAxisForCollision(Eigen::Vector3d position, double distance_to_be_inside, int axis){

	std::vector<int> spheres_id_to_check;
	for (auto i = 0; i < spheres.size(); ++i) {

            double min_i = spheres[i].P[axis] - spheres[i].radius;

			if (min_i> position[axis] + distance_to_be_inside) {
				continue;
			}
			else {
				double max_i = spheres[i].P[axis] + spheres[i].radius;

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
std::vector<int> findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3) {
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

bool Axon::isPosInsideAxon_(Eigen::Vector3d position, double distance_to_be_inside){
  
    //cout << "isSphereInsideAxon_ : " << id << endl;
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

                return true;
            }
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

    Eigen::Vector2d x_limits = {boxes[i].x_min - distance_to_be_inside, boxes[i].x_max + distance_to_be_inside};
    Eigen::Vector2d y_limits = {boxes[i].y_min - distance_to_be_inside, boxes[i].y_max + distance_to_be_inside};
    Eigen::Vector2d z_limits = {boxes[i].z_min - distance_to_be_inside, boxes[i].z_max + distance_to_be_inside};
    for (int j = 0; j < 3; ++j) {
        if (position[j] < x_limits[0] || position[j] > x_limits[1] ) {
            return false;
        }
        if (position[j] < y_limits[0] || position[j] > y_limits[1]) {
            return false;
        }
        if (position[j] < z_limits[0] || position[j] > z_limits[1] ) {
            return false;
        }
    }
    return true;
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