#include "Axon.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>


using namespace Eigen;
using namespace std;


Axon::Axon(const Axon &ax)
{

    id = ax.id;
    spheres = ax.spheres;
    radius = ax.radius;
    begin = ax.begin;
    end = ax.end;
    Box = ax.Box;
    //projections = ax.projections;

}



void Axon::set_spheres(std::vector<Sphere> spheres_to_add){

    // value of center of sphere at x that has the highest x center value
    double sph_highest_x_val = spheres_to_add[0].P[0] + spheres_to_add[0].radius;
    // value of center of sphere at y that has the highest y center value
    double sph_highest_y_val= spheres_to_add[0].P[1] +spheres_to_add[0].radius;
    // value of center of sphere at x that has the lowest x center value
    double sph_lowest_x_val= spheres_to_add[0].P[0]- spheres_to_add[0].radius;
    // value of center of sphere at y has the lowest y center value
    double sph_lowest_y_val= spheres_to_add[0].P[1]- spheres_to_add[0].radius;



    // set axon_id of spheres to the id of axon
    for (int i=0; i<spheres_to_add.size(); ++i){
        spheres_to_add[i].ax_id = id;
        spheres_to_add[i].id = i;

        // highest y value
        if (spheres_to_add[i].P[1]+spheres_to_add[i].radius > sph_highest_y_val){
     
            sph_highest_y_val = spheres_to_add[i].P[1]+spheres_to_add[i].radius;
        }
        // highest x value
        if (spheres_to_add[i].P[0]+spheres_to_add[i].radius > sph_highest_x_val){
 
            sph_highest_x_val = spheres_to_add[i].P[0]+spheres_to_add[i].radius;
        }
        // lowest x value
        if (spheres_to_add[i].P[0]-spheres_to_add[i].radius < sph_lowest_x_val){
     
            sph_lowest_x_val = spheres_to_add[i].P[0]-spheres_to_add[i].radius;
        }
        // lowest y value
        if (spheres_to_add[i].P[1]-spheres_to_add[i].radius  < sph_lowest_y_val){
       
            sph_lowest_y_val = spheres_to_add[i].P[1]-spheres_to_add[i].radius;
        }
    
        //cout << "projections.sph_projections_x.size() in set_spheres: "<< projections.sph_projections_x.size() << endl;
        

    }
    if (spheres_to_add.size() != 0){

        this->begin = spheres_to_add[0].P;

        this->spheres = spheres_to_add;

        this->end = spheres_to_add[spheres_to_add.size()-1].P;

        Box.clear();
        //x
        Box.push_back({sph_lowest_x_val , sph_highest_x_val});
        //y
        Box.push_back({sph_lowest_y_val , sph_highest_y_val});
        //z
        Box.push_back({spheres_to_add[0].P[2] - spheres_to_add[0].radius, spheres_to_add[spheres_to_add.size()-1].P[2] + spheres_to_add[spheres_to_add.size()-1].radius});
    }
    // add projection of spheres 
    //for (int i=0; i<spheres_to_add.size(); ++i){
    //    add_projection(spheres_to_add[i]);
    //}
}
/*
void Axon::add_projection(Sphere sphere_to_add){

    // 2d limits of axon
    Vector2d x_limits;
    Vector2d y_limits;
    // projections are in descending order. When added, it is added at the right position.
    for (unsigned axis = 0; axis < 3; ++axis) {

        double position1, position2;

        // projections
        // P + radius
        position1 = sphere_to_add.P[axis] + sphere_to_add.radius;
        Projections::projection_pt p1 {position1, id, sphere_to_add.id};
        // P - radius
        position2 = sphere_to_add.P[axis] - sphere_to_add.radius;
        Projections::projection_pt p2 {position2, id, sphere_to_add.id};

        projections.append_right_place(p1, p2, axis);
        //cout << "projections.sph_projections_x.size() in add_proj: "<< projections.sph_projections_x.size() << endl;
        
    }

}
*/

bool check_with_edge(Eigen::Vector3d position, Eigen::Vector2d x_limits, Eigen::Vector2d y_limits){ 
    if ((position[0] >=  x_limits[0])  && (position[0] <= x_limits[1])){
        if ((position[1] >= y_limits[0]) && (position[1] <=  y_limits[1])){
            return true;
        }
    }
    return false;
} 



bool Axon::intersection_sphere_vector(double &t1, double &t2, Sphere &s, Eigen::Vector3d &step, double &step_length, Eigen::Vector3d &pos, double &c){
    //https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    Eigen::Vector3d m = pos - s.P;
    double rad = s.radius;

    double a = 1;
    double b = (m.dot(step));
    c = m.dot(m) - rad*rad;
    double discr = b*b - a*c;

    if (discr < 0.0 ){
        return false;
    }

    t1 = (-b + sqrt(discr))/(a);
    t2 = (-b - sqrt(discr))/(a);

    return true;

}


bool check_inside(vector<double> list_c){
    for (unsigned i=0 ; i< list_c.size(); ++i){
        if (list_c[i] < 1e-10){
            return true;
            break;        
        }
    }
    return false;
}

bool check_outside(vector<double> list_c){
    for (unsigned i=0 ; i< list_c.size(); ++i){
        if (list_c[i] < -1e-10){
            return false;      
        }
    }
    return true;
}


bool Axon::checkCollision(Walker &walker,  Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{
    
    string message;
    Eigen::Vector3d O;
    walker.getVoxelPosition(O);

    // distances to intersections
    std::vector<double> dist_intersections;
    // values indicating whether the walker is inside or outside a sphere
    std::vector<double> cs;
    std::vector<int> sph_ids;


    int first_ind = 0;
    int last_ind = spheres.size();

    //cout << "walker.in_sph_index.size() :"<< walker.in_sph_index.size()<< endl;
    /*
    if (walker.in_sph_index.size()>0){
        first_ind = walker.in_sph_index[0];
        last_ind = walker.in_sph_index[walker.in_sph_index.size()-1];
    }
    else{
        first_ind = last_ind = 0;
    }
    */

    for (unsigned i= first_ind ; i< last_ind; ++i){
        
        // distances to collision
        double t1;
        double t2;
        double c;
        bool intersect = intersection_sphere_vector(t1, t2, spheres[i], step, step_lenght, O, c); 

        if (intersect){


            Eigen::Vector3d position; 
            std::vector<int> sph_ids_;
            bool isinside; 
            
            //if the collision are too close or negative.
            if(Walker::bouncing){
                if( t1 >= EPS_VAL ){
                    
                    if (walker.initial_location== Walker::intra){
                        position = walker.pos_v + t1*step;
                        isinside = isPosInsideAxon_(position, -barrier_tickness,sph_ids_);
                        if (!isinside || sph_ids_.size()<=1){
                            dist_intersections.push_back(t1);
                            sph_ids.push_back(i);
                            cs.push_back(c);
                        }
                    }
                    else{

                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                        
                    }
                }
                
                if(t2 >= EPS_VAL ){
     
                    if (walker.initial_location== Walker::intra){
                        position = walker.pos_v + t2*step;
                        isinside = isPosInsideAxon_(position, -barrier_tickness,sph_ids_);
                        if (!isinside || sph_ids_.size()<=1){
                            dist_intersections.push_back(t2);
                            sph_ids.push_back(i);
                            cs.push_back(c);
                        }
                    }
                    else{

                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                        
                    }
                }
                
            }
            else{
                if( t1 >= 0 ){
                    
                    if (walker.initial_location== Walker::intra){
                        position = walker.pos_v + t1*step;
                        isinside = isPosInsideAxon_(position, -barrier_tickness,sph_ids_);
                        if (!isinside || sph_ids_.size()<=1){
                            dist_intersections.push_back(t1);
                            sph_ids.push_back(i);
                            cs.push_back(c);
                        }
                    }
                    else{

                        dist_intersections.push_back(t1);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                        
                    }
                }
                
                if(t2 >= 0){
                    
                    if (walker.initial_location== Walker::intra){
                        position = walker.pos_v + t2*step;
                        isinside = isPosInsideAxon_(position, -barrier_tickness,sph_ids_);
                        if (!isinside || sph_ids_.size()<=1){
                            dist_intersections.push_back(t2);
                            sph_ids.push_back(i);
                            cs.push_back(c);
                        }
                    }
                    else{

                        dist_intersections.push_back(t2);
                        sph_ids.push_back(i);
                        cs.push_back(c);
                        
                    }
                    
                }  
                
            }
            
        } 
    }

    if(dist_intersections.size() > 0){


        unsigned index_;

        auto min_distance_int = std::min_element(std::begin(dist_intersections), std::end(dist_intersections));
        index_ = std::distance(std::begin(dist_intersections), min_distance_int);
        

        int sphere_ind = sph_ids[index_];

        double dist_to_collision = dist_intersections[index_];


        if (dist_to_collision <= step_lenght + barrier_tickness){

            //Eigen::Vector3d position = walker.pos_v + dist_to_collision*step;
            //std::vector<int> sph_ids_;

            //bool isinside = isPosInsideAxon(position,  barrier_tickness, radius+barrier_tickness, sph_ids_);
                    
            //cout << "number of spheres position is inside : " << sph_ids_.size() << endl;

            // if extra, cannot hit sphere from inside
            if (walker.initial_location== Walker::extra){
                if (cs[index_]< -1e-15){
                    colision.col_location = Collision::inside;
                    walker.location = Walker::intra;
                    //cout << cs[index_] << endl;


                    //Sphere s = spheres[sphere_ind];
                    //cout << "distance to walker : " << s.minDistance(O)<<endl;
                    //cout << "distance to collision point : " << s.minDistance(O+step*dist_to_collision)<<endl;
                }
                else if (cs[index_]> 1e-15){
                    colision.col_location = Collision::outside;
                }
                else{
                    colision.col_location = Collision::unknown;
                    //cout << "UNKNOWN" << endl;
                    //cout << cs[index_] << endl;

                    //Sphere s = spheres[sphere_ind];
                    //cout << "distance to walker : " << s.minDistance(O)<<endl;
                    //cout << "distance to collision point : " << s.minDistance(O+step*dist_to_collision)<<endl;
                
                }
            }
            else{
                colision.col_location = Collision::inside;
            }

            colision.type = Collision::hit;
            colision.rn = cs[index_];  
            
    
            //std::cout << "Collision, sphere index :" << sph_ids[index_] << endl;

            //std::cout << "Collision, walker position :" << O << endl;
            //std::cout << "Collision, collision position :" << O+step*dist_to_collision << endl;
            
            colision.t = fmin(dist_to_collision,step_lenght);
            colision.colision_point = walker.pos_v + colision.t*step;
            Eigen::Vector3d normal = (colision.colision_point- spheres[sphere_ind].P).normalized();
            
            Eigen::Vector3d temp_step = step;

            elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);
            
            colision.bounced_direction = temp_step.normalized();

            return true;
        }
        else{
            colision.type = Collision::null;
            return false;   
        }
        
    }
    else{
        /*
        Eigen::Vector3d position = walker.pos_v + step_lenght*step;
        std::vector<int> sph_ids_;
        bool isinside = isPosInsideAxon_(position, -2*barrier_tickness, sph_ids_);
        if (isinside){
            cout << " position after not bouncing is inside" << endl;
        }
        */
        
        colision.type = Collision::null;
        return false;   
    }

}

bool Axon::isWalkerInsideAxon(Walker &walker, double distance_to_be_inside){
    
    Eigen::Vector3d O;
    walker.getVoxelPosition(O);
    std::vector<int> sph_ids;
    bool isinside = isPosInsideAxon_(O, distance_to_be_inside, sph_ids);
    
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

bool Axon::isPosInsideAxon_(Eigen::Vector3d position, double distance_to_be_inside, std::vector<int> &sph_ids){
    sph_ids.clear();
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
                sph_ids.push_back(spheres_to_check_all_axes[i]);

                //return true;
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
    
    //return false;
}
/*
bool Axon::isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside, double max_radius, std::vector<int> &sph_ids){
    // when checking collision with walker -> check with normal radius
    // when checking with collisions of other axons -> check with max_radius so there is room for swelling
    std::vector<std::vector<Projections::projection_pt>> coliding_projs;
    bool colliding_all_axes;
    Sphere sphere_ ;
    double rad;
    sph_ids.clear();

    // if position is in box with axon inside
    if(isNearAxon(position, distance_to_be_inside)){
 
        // find all projections in between the two projections of the edges
        //cout << "is near " << endl;
        coliding_projs = projections.find_collisions_all_axes(position, max_radius + barrier_tickness, id, distance_to_be_inside);

        if (coliding_projs.size() == 3){ 
            
            // for all coliding objects in x 
            //cout << "debug " << coliding_projs[0].size() << endl;

            for(unsigned j = 0; j < coliding_projs[0].size() ; j++){ 

                const Projections::projection_pt coliding_proj = coliding_projs[0][j];
                // if the same coliding objects are also in y and z but are not from same objects

                colliding_all_axes = (projections.isProjInside(coliding_projs[1], coliding_proj) && projections.isProjInside(coliding_projs[2], coliding_proj));

                if (colliding_all_axes){
          
                    sphere_ = spheres[coliding_proj.sph_id];
                    //cout << "coliding_projs[0][j].sph_id :" << coliding_projs[0][j].sph_id << endl;
              
                    if (sphere_.minDistance(position) <= distance_to_be_inside){ 
                        //cout << "is near sphere "<< coliding_proj.sph_id << endl;
                        sph_ids.push_back(sphere_.id);
                    
                    }  
                }
            }
        }
        
    }
    if (sph_ids.size()>0){
        // Sort sph_ids in ascending order
        std::sort(sph_ids.begin(), sph_ids.end());

        return true;
    }
    else{
        return false;
    }
} 
*/

bool Axon::isNearAxon(Eigen::Vector3d position, double distance_to_be_inside){
    bool isnear = false;
    Eigen::Vector2d x_limits = Box[0];
    Eigen::Vector2d y_limits = Box[1];

    if(check_with_edge(position, x_limits+ Eigen::Vector2d{-distance_to_be_inside,distance_to_be_inside}, y_limits+ Eigen::Vector2d{-distance_to_be_inside,distance_to_be_inside})){
       return true;
    }  
    /*
    // if the axon touches the border of the voxel, we must take into account the translation
    // of the axon at the opposite plane
    if (x_limits[0]<0){
        if(check_with_edge(position,{x_limits[0]+ end[2] - distance_to_be_inside, x_limits[1] + end[2]+distance_to_be_inside} , y_limits)){
            return true;
        } 
    }
    else if (x_limits[0]>end[2]){
        if(check_with_edge(position,{x_limits[0] - end[2]- distance_to_be_inside, x_limits[1] - end[2]+distance_to_be_inside} , y_limits)){
            return true;
        } 
    } 

    if (y_limits[0]<0){
        if(check_with_edge(position,x_limits,{y_limits[0]+ end[2]- distance_to_be_inside, y_limits[1] + end[2]+distance_to_be_inside} )){
            return true;
        } 
    }
    else if (y_limits[0]>end[2]){
        if(check_with_edge(position,x_limits, {y_limits[0] - end[2]- distance_to_be_inside, y_limits[1] - end[2]+distance_to_be_inside})){
            return true;
        } 
    } 
    */

    return false;
}

bool Axon::isNearAxon(Walker walker, double distance_to_be_inside){
    bool isnear = false;
    Eigen::Vector2d x_limits = Box[0];
    Eigen::Vector2d y_limits = Box[1];

    Eigen::Vector3d position;
    walker.getVoxelPosition(position);

    return isNearAxon(position,distance_to_be_inside);

}



