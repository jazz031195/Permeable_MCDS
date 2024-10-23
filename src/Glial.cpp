#include "Glial.h"
#include "Eigen/Dense"
#include "constants.h"
#include <numeric>

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
    Box = gl.Box;

};

void Glial::set_spheres(const std::vector<Sphere> &spheres_to_add){

    for (int i = 0; i < spheres_to_add.size(); i++){

        Sphere sphere_to_add = spheres_to_add[i];
        // value of center of sphere at x that has the highest x center value
        double sph_highest_x_val;
        // value of center of sphere at y that has the highest y center value
        double sph_highest_y_val;
        // value of center of sphere at x that has the lowest x center value
        double sph_lowest_x_val;
        // value of center of sphere at y has the lowest y center value
        double sph_lowest_y_val;
        // value of center of sphere at z that has the highest z center value
        double sph_highest_z_val;
        // value of center of sphere at z that has the lowest z center value
        double sph_lowest_z_val;

        if (this->Box.empty() && soma.radius > 0){
            // intialise box to soma
            // create box around that one sphere
            sph_highest_x_val = soma.P[0]+ 2.0*soma.radius;
            sph_lowest_x_val = soma.P[0] -2.0*soma.radius;
            sph_highest_y_val = soma.P[1] +2.0*soma.radius;
            sph_lowest_y_val = soma.P[1] -2.0*soma.radius;
            sph_highest_z_val = soma.P[2] +2.0*soma.radius;
            sph_lowest_z_val = soma.P[2] -2.0*soma.radius;
            //x
            Box.push_back({sph_lowest_x_val, sph_highest_x_val});
            //y
            Box.push_back({sph_lowest_y_val, sph_highest_y_val});
            //z
            Box.push_back({sph_lowest_z_val, sph_highest_z_val});

        }
        
        else{
            // take values in the box
            sph_highest_x_val = Box[0][1];
            sph_lowest_x_val = Box[0][0];
            sph_highest_y_val = Box[1][1];
            sph_lowest_y_val = Box[1][0];
            sph_highest_z_val = Box[2][1];
            sph_lowest_z_val = Box[2][0];

        }
        //cout << "---------" << endl;
        //cout << "Box: " << Box[0][0] << " " << Box[0][1] << " " << Box[1][0] << " " << Box[1][1] << " " << Box[2][0] << " " << Box[2][1] << endl;

        // if the extemity of sphere is higher than extermity of Box (x)
        if (sph_highest_x_val < sphere_to_add.P[0]+ 2.0*sphere_to_add.radius){
            Box[0][1] = sphere_to_add.P[0]+ 2.0*sphere_to_add.radius;
        }
        // if the extemity of sphere is lower than extermity of Box (x)
        if(sph_lowest_x_val > sphere_to_add.P[0]- 2.0*sphere_to_add.radius){
            Box[0][0] = sphere_to_add.P[0]- 2.0*sphere_to_add.radius;
        }
        // if the extemity of sphere is higher than extermity of Box (y)
        if (sph_highest_y_val < sphere_to_add.P[1]+ 2.0*sphere_to_add.radius){
            Box[1][1] = sphere_to_add.P[1]+ 2.0*sphere_to_add.radius;
        }
        // if the extemity of sphere is lower than extermity of Box (y)
        if(sph_lowest_y_val > sphere_to_add.P[1]- sphere_to_add.radius){
            Box[1][0]= sphere_to_add.P[1]- 2.0*sphere_to_add.radius;
        }
        // if the extemity of sphere is higher than extermity of Box (z)
        if (sph_highest_z_val < sphere_to_add.P[2]+ 2.0*sphere_to_add.radius){
            Box[2][1] = sphere_to_add.P[2]+ 2.0*sphere_to_add.radius;
        }
        // if the extemity of sphere is lower than extermity of Box (z)
        if(sph_lowest_z_val > sphere_to_add.P[2]- 2.0*sphere_to_add.radius){
            Box[2][0] = sphere_to_add.P[2]- 2.0*sphere_to_add.radius;
        }
    //cout << "Box: " << Box[0][0] << " " << Box[0][1] << " " << Box[1][0] << " " << Box[1][1] << " " << Box[2][0] << " " << Box[2][1] << endl;

    }
    processes = spheres_to_add;

    
}


bool Glial::isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside){
    
    if (Box.empty()){
        return false;
    }
    if ((position[0] >=  Box[0][0]- distance_to_be_inside)  && (position[0] <= Box[0][1] + distance_to_be_inside)){
        if ((position[1] >= Box[1][0] - distance_to_be_inside) && (position[1] <=  Box[1][1] + distance_to_be_inside)){
            return true;
        }
    }
    return false;
}

std::vector<int> Glial::checkAxisForCollision(Eigen::Vector3d position, double distance_to_be_inside, int axis){

	std::vector<int> spheres_id_to_check;
    std::vector<Sphere> all_spheres;
    if (soma.radius > 0)
        all_spheres = {soma};
    
    all_spheres.insert(all_spheres.end(), processes.begin(), processes.end());

    //cout << "all_spheres.size() : " << all_spheres.size() << endl;
	for (auto i = 0; i < all_spheres.size(); ++i) {

            double min_i = all_spheres[i].P[axis] - all_spheres[i].radius;
			if (min_i - distance_to_be_inside > position[axis]) {
				continue;
			}
			else {
				double max_i = all_spheres[i].P[axis] + all_spheres[i].radius;
                if (position[axis] > max_i + distance_to_be_inside) {
                    continue;
                }
                else{
                    if ((soma.radius > 0) && (i > 0))
                        spheres_id_to_check.push_back(i-1);
                    else
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

bool Glial::isPosInsideGlialCell_(const Eigen::Vector3d& position, const double& distance_to_be_inside){
    for (auto i = 0; i < processes.size(); ++i) {
        if (processes[i].minDistance(position) <= distance_to_be_inside){
            return true;
        }
        //cout <<"processes[i].minDistance(position) : " << processes[i].minDistance(position)  << endl;
    }
    if (soma.minDistance(position) <= distance_to_be_inside){
        
        return true;
    }
    //cout <<"soma.minDistance(position) : " << soma.minDistance(position) << endl;
    return false;
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
                Sphere sphere_to_check = processes[spheres_to_check_all_axes[i]];
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
    bool isinside = FindSphereinGlial(pos, step_lenght,sph_ids_walker_is_inside);

    for (auto i = 0; i < sph_ids_walker_is_inside.size() ; ++i) {


        Sphere sphere_to_check;

        if (sph_ids_walker_is_inside[i] == 0){
            //soma
            sphere_to_check = soma;
        }
        else{
            //processes
            sphere_to_check = processes[sph_ids_walker_is_inside[i]];
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

bool Glial::checkCollision(Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, Collision &colision){
    // distances to intersections
    std::vector<double> dist_intersections;

    //all_spheres.insert(all_spheres.end(), processes.begin(), processes.end());

    walker.previous_location = walker.location;

    //cout << "-------------------" << endl; ;

    std::vector<int> spheres_ids;

    //cout << "all_spheres.size() : " << all_spheres.size() << endl;

    //cout << "walker.sph_id_to_check size : " << walker.sph_id_to_check.size() << endl;

    find_all_intersections(walker, step, step_lenght + barrier_tickness, dist_intersections, spheres_ids);
    
    if (dist_intersections.empty()){
        //cout << "dist_intersections.empty()" << endl;
        if (walker.location == Walker::intra){
            bool isinside = isPosInsideGlialCell(walker.pos_v, EPS_VAL);
            if (!isinside){
                colision.col_location = Collision::outside;
                walker.in_obj_index = -1;
                walker.in_obj_type = -1;
                walker.location = Walker::extra;
                cout << "no intersection but outside " << endl;
                return true;
            }
        }
        colision.type = Collision::null;
        return false;

    }
    else{

        // sort indexes based on dist_intersections
        std::vector<int> idx(dist_intersections.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&dist_intersections](int i1, int i2) {return dist_intersections[i1] < dist_intersections[i2];});


        //bool isnearEdge_walker = !isPosInsideGlialCell_(walker.pos_v , -EPS_VAL);
        //for (auto i = 0; i < idx.size(); ++i) {
        //    int index = idx[i];
        //    Eigen::Vector3d pos = walker.pos_v + dist_intersections[index]*step;
        //    bool isnearEdge_minus = !isPosInsideGlialCell_(pos, -EPS_VAL);
        //    cout << " step length :" <<step_lenght <<" i : " <<i <<" index : " << index << " dist_intersections[i] : " << dist_intersections[index] << " id_to_check[i] : " <<  spheres_ids[index] << " isnearEdge_walker : " << isnearEdge_walker << " isnearEdge_minus : " << isnearEdge_minus  << endl;
        //}

        for (auto i = 0; i < idx.size(); ++i) {
            int index = idx[i];
            Eigen::Vector3d pos = walker.pos_v + dist_intersections[index]*step;
            bool isnearEdge = true;
            if (walker.location == Walker::intra){  
                isnearEdge = !isPosInsideGlialCell(pos, -EPS_VAL);
            }

            //bool isnearEdge = true;

            if (isnearEdge ){
                //cout <<"Chosen: step length :" <<step_lenght <<" i : " <<i <<" index : " << index << " dist_intersections[i] : " << dist_intersections[index] << " id_to_check[i] : " <<  spheres_ids[index] << " isnearEdge : " << isnearEdge << endl;

                if (dist_intersections[index] <= step_lenght+barrier_tickness ){
                    

                    Sphere sph;
                    if (spheres_ids[index] == 0){
                        sph = soma;
                    }
                    else{
                        sph = processes[spheres_ids[index]];
                    }
                                
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
                    //cout << "bounced_step : " << step << " step : "<< temp_step << endl;
                    //cout << " collision point : " << colision.colision_point <<endl;

                    colision.bounced_direction = temp_step.normalized();
                    colision.perm_crossing = 0.;

                    // if(walker.location == Walker::extra)
                    // {

                    // }
                    if (rn < -1e-10){
                        colision.col_location = Collision::inside;
                        walker.in_obj_index = id;
                        walker.in_obj_type = 1;
                        walker.location = Walker::intra;
                        
                    }
                    else if (rn > 1e-10){
                        colision.col_location = Collision::outside;
                        walker.in_obj_index = -1;
                        walker.in_obj_type = -1;
                        walker.location = Walker::extra;
                 
                        
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

                    //bool isnearEdge_walker = !isPosInsideGlialCell(walker.pos_v, -barrier_tickness);

                    //if (isnearEdge_walker){
                    //    colision.col_location = Collision::on_edge;
                    //    colision.bounced_direction = -step;
                    //    cout <<"ON EDGE" << endl;
                    //}

                    colision.perm_crossing = 0.;
                    //cout <<"hit" << endl;
            
                    return true;
                }
                else{
                    //cout << "too far"<< endl;
                    colision.type = Collision::null;
                 
                    return false;
                }
            }
            
        }
        colision.type = Collision::null;
      
        return false;
        
    }

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
                Sphere sphere_to_check = processes[spheres_to_check_all_axes[i]];
                if (sphere_to_check.minDistance(position) <= distance_to_be_inside){
                    sph_ids.push_back(spheres_to_check_all_axes[i]);
                }
            }
        }
        spheres_id_to_check.clear();
        spheres_to_check_all_axes.clear();
    }


    if (sph_ids.size()>0){
        // Sort sph_ids in ascending order
        //std::sort(sph_ids.begin(), sph_ids.end());
        return true;
    }
    else{
        return false;
    }
    
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


double Glial::minDistance(Walker &w){
    //Origin of the ray
    Vector3d O;
    w.getVoxelPosition(O);
    // find distance to Box
    double minDistSquared = 0.0;

    for (int i = 0; i < 3; ++i) {
        double v = O[i];
        double min = Box[i][0];
        double max = Box[i][1];

        if (v < min) {
            minDistSquared += (min - v) * (min - v);
        } else if (v > max) {
            minDistSquared += (v - max) * (v - max);
        }
    }

    return std::sqrt(minDistSquared);
}
