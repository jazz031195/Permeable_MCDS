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
    Box_branch = gl.Box_branch;

};

void Glial::init_box(){

    double sph_highest_x_val = soma.P[0]+ soma.radius;
    double sph_lowest_x_val  = soma.P[0] -soma.radius;
    double sph_highest_y_val = soma.P[1] +soma.radius;
    double sph_lowest_y_val  = soma.P[1] -soma.radius;
    double sph_highest_z_val = soma.P[2] +soma.radius;
    double sph_lowest_z_val  = soma.P[2] -soma.radius;
    //x
    Box.push_back({sph_lowest_x_val, sph_highest_x_val});
    //y
    Box.push_back({sph_lowest_y_val, sph_highest_y_val});
    //z
    Box.push_back({sph_lowest_z_val, sph_highest_z_val});
}

void Glial::set_box_neuron(){

    init_box();
    for(int axis = 0; axis < 3; axis++)
    {
        double sph_highest = soma.P[axis] + 2.0*soma.radius;
        double sph_lowest  = soma.P[axis] - 2.0*soma.radius;

        for (size_t i = 0; i < Box_branch.size(); i++)
        {
            auto box = Box_branch[i][axis];
            if (box[1] > sph_highest)
                sph_highest = box[1];
            // if the extemity of sphere is lower than extermity of Box (x)
            if(box[0] < sph_lowest)
                sph_lowest = box[0];
        } 
        Box[axis] = {sph_lowest, sph_highest};  
    } 
}

void Glial::set_spheres(const std::vector<Sphere> &spheres_to_add){

    std::vector<Eigen::Vector2d> Box_branch_all_axis;
    for (int axis = 0; axis < 3; axis++)
    {
        for (int i = 0; i < spheres_to_add.size(); i++){

            Sphere sphere_to_add = spheres_to_add[i];

            double sph_highest, sph_lowest;
            sph_highest = sphere_to_add.P[axis] + 2.0*sphere_to_add.radius;
            sph_lowest  = sphere_to_add.P[axis] - 2.0*sphere_to_add.radius;

            if (i == 0)
            {
                // intialise box to soma
                // create box around that one sphere
                Box_branch_all_axis.push_back({sph_lowest, sph_highest});
            }
            else
            {
                // if the extemity of sphere is higher than extermity of Box (x)
                if (sph_highest > Box_branch_all_axis[axis][1])
                    Box_branch_all_axis[axis][1] = sph_highest;
                if (sph_lowest < Box_branch_all_axis[axis][0])
                    Box_branch_all_axis[axis][0] = sph_lowest;

            }
        }      
    }

    Box_branch.push_back(Box_branch_all_axis);
    processes.push_back(spheres_to_add);
}


bool Glial::isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside){
    
    if (Box.empty())
        return false;
    
    for(int i=0; i < 3; i++)
    {   
        if ((position[i] >=  Box[i][0] - distance_to_be_inside)  && (position[i] <= Box[i][1] + distance_to_be_inside))
            continue;
        else
            return false;
    }

    return true;
}

bool Glial::isNearDendrite(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int>& dendrite_ids){
    
    if (Box_branch.empty())
        return false;
    
    
    for(int i=0; i < Box_branch.size(); i++)
    {
        int count = 0;
        for(int axis=0; axis < 3; axis++)
        {
            if ((position[axis] >=  Box_branch[i][axis][0] - distance_to_be_inside)  && (position[axis] <= Box_branch[i][axis][1] + distance_to_be_inside))
                count++;
        }
        if (count == 3)
            dendrite_ids.push_back(i);
    }

    return dendrite_ids.size() > 0;
}

std::vector<int> Glial::checkAxisForCollision_soma(Eigen::Vector3d position, double distance_to_be_inside, int axis){

	std::vector<int> spheres_id_to_check;
    if (soma.radius > 0)
    {
        double min_i = soma.P[axis] - soma.radius;
        double max_i = soma.P[axis] + soma.radius;
        // cout << "soma " << min_i - distance_to_be_inside << " < " << position[axis] << endl;
        // cout << "soma " << position[axis] << " <= " << max_i + distance_to_be_inside << endl;
        if ((min_i - distance_to_be_inside < position[axis]) && (position[axis] <= max_i + distance_to_be_inside))
            spheres_id_to_check.push_back(0);
    }
    return spheres_id_to_check;
}

std::vector<int> Glial::checkAxisForCollision_dendrite(Eigen::Vector3d position, double distance_to_be_inside, int axis, int dendrite_to_check){

	std::vector<int> spheres_id_to_check;
    std::vector<Sphere> all_spheres;


    all_spheres.insert(all_spheres.end(), processes[dendrite_to_check].begin(), processes[dendrite_to_check].end());

	for (auto i = 0; i < all_spheres.size(); ++i) {

        double min_i = all_spheres[i].P[axis] - all_spheres[i].radius;
        double max_i = all_spheres[i].P[axis] + all_spheres[i].radius;
        // cout << "de " << min_i - distance_to_be_inside << " < " << position[axis] << endl;
        // cout << "de " << position[axis] << " <= " << max_i + distance_to_be_inside << endl;
        if ((min_i - distance_to_be_inside < position[axis]) && (position[axis] <= max_i + distance_to_be_inside))
            spheres_id_to_check.push_back(i);
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

    if (isNearGlialCell(position, distance_to_be_inside))
    {
        std::vector<std::vector<int>> spheres_id_to_check;
        
        // check soma
        for (auto axis = 0; axis < 3; ++axis) 
        {
            int count = 0;
            for (auto axis = 0; axis < 3; ++axis)
            {
                spheres_id_to_check.push_back(checkAxisForCollision_soma(position, distance_to_be_inside, axis)); // check for collision along 1 axis
                if(spheres_id_to_check[spheres_id_to_check.size() - 1].size() > 0)
                    ++count;
            }

            if(count == 3)
            {

                // find common ids in all 3 axes
                std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
                for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) 
                {
                    if (soma.minDistance(position) <= distance_to_be_inside)
                    {
                        return true;
                    }
                }
                spheres_to_check_all_axes.clear();
            }
        }

        spheres_id_to_check.clear();
        // check dendrites
        std::vector<int> dendrite_ids;
        bool isNearDendrite_ = isNearDendrite(position, distance_to_be_inside, dendrite_ids);
        if(isNearDendrite_)
        {

            for (auto d_id = 0; d_id < dendrite_ids.size(); ++d_id) 
            {
                int count = 0;
                for (auto axis = 0; axis < 3; ++axis)
                {
                    spheres_id_to_check.push_back(checkAxisForCollision_dendrite(position,distance_to_be_inside, axis, dendrite_ids[d_id])); // check for collision along 1 axis
                    if(spheres_id_to_check[spheres_id_to_check.size() - 1].size() > 0)
                    ++count;
                }
                    
                if(count == 3)
                {
                    // find common ids in all 3 axes
                    std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
                    for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) 
                    {
                        Sphere sphere_to_check = processes[dendrite_ids[d_id]][spheres_to_check_all_axes[i]];
                        if (sphere_to_check.minDistance(position) <= distance_to_be_inside)
                        {
                            return true;

                        }
                    }
                    spheres_to_check_all_axes.clear();
                }
                spheres_id_to_check.clear();
            }
        }
        return false;
    }
    else
        return false;
}

bool Glial::isPosInsideGlialCell_verbose(const Eigen::Vector3d& position, const double& distance_to_be_inside){
    cout << "is pos inside glial? " << endl;

    if (isNearGlialCell(position, distance_to_be_inside))
    {
        cout << "is near glial " << endl;
        std::vector<std::vector<int>> spheres_id_to_check;
        
        // check soma
        for (auto axis = 0; axis < 3; ++axis) 
        {
            int count = 0;
            for (auto axis = 0; axis < 3; ++axis)
            {
                spheres_id_to_check.push_back(checkAxisForCollision_soma(position, distance_to_be_inside, axis)); // check for collision along 1 axis
                if(spheres_id_to_check[spheres_id_to_check.size() - 1].size() > 0)
                    ++count;
            }

            if(count == 3)
            {
                cout << "is near soma " << endl;

                // find common ids in all 3 axes
                std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
                for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) 
                {
                    if (soma.minDistance(position) <= distance_to_be_inside)
                    {
                        cout << "is in soma " << endl;
                        return true;
                    }
                }
                spheres_to_check_all_axes.clear();
            }
        }

        spheres_id_to_check.clear();
        // check dendrites
        std::vector<int> dendrite_ids;
        bool isNearDendrite_ = isNearDendrite(position, distance_to_be_inside, dendrite_ids);
        if(isNearDendrite_)
        {
            cout << "is near dendrite " << endl;

            for (auto d_id = 0; d_id < dendrite_ids.size(); ++d_id) 
            {
                int count = 0;
                for (auto axis = 0; axis < 3; ++axis)
                {
                    spheres_id_to_check.push_back(checkAxisForCollision_dendrite(position,distance_to_be_inside, axis, dendrite_ids[d_id])); // check for collision along 1 axis
                    if(spheres_id_to_check[spheres_id_to_check.size() - 1].size() > 0)
                    ++count;
                }
                    
                if(count == 3)
                {
                    cout << "is near dendrite 2" << endl;

                    // find common ids in all 3 axes
                    std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
                    for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) 
                    {
                        Sphere sphere_to_check = processes[dendrite_ids[d_id]][spheres_to_check_all_axes[i]];
                        cout << sphere_to_check.minDistance(position) << " <= " << distance_to_be_inside << endl;
                        if (sphere_to_check.minDistance(position) <= distance_to_be_inside)
                        {
                            cout << "is in dendrite " << dendrite_ids[d_id] << " sph " << spheres_to_check_all_axes[i] << endl;
                            return true;

                        }
                    }
                    spheres_to_check_all_axes.clear();
                }
                spheres_id_to_check.clear();
            }
        }
        return false;
    }
    else
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

void Glial::find_all_intersections(const Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, std::vector<double>& dist_intersections, std::vector<Sphere*>& spheres){

    dist_intersections.clear();
    spheres.clear();

    std::vector<Sphere*> sph_walker_is_inside;
    Eigen::Vector3d pos = walker.pos_v;
    bool isinside = FindSphereinGlial(pos, step_lenght, sph_walker_is_inside);

    for (auto i = 0; i < sph_walker_is_inside.size() ; ++i) 
    {

        if (sph_walker_is_inside[i]->minDistance(walker.pos_v) > step_lenght)
            continue;

        double t1, t2;
        bool intersect = intersection_sphere_vector(t1, t2, *sph_walker_is_inside[i], step, step_lenght, walker.pos_v);
        if (intersect)
        {
            if (walker.status == Walker::bouncing){
                if (t1 > EPS_VAL)
                {
                    dist_intersections.push_back(t1);
                    spheres.push_back(sph_walker_is_inside[i]);
                }
                if (t2 > EPS_VAL)
                {
                    dist_intersections.push_back(t2);
                    spheres.push_back(sph_walker_is_inside[i]);
                }
            }
            else{
                if (t1 > 0)
                {
                    dist_intersections.push_back(t1);
                    spheres.push_back(sph_walker_is_inside[i]);
                }
                if (t2 > 0)
                {
                    dist_intersections.push_back(t2);
                    spheres.push_back(sph_walker_is_inside[i]);
                }
            }
        }
    }
}

bool Glial::checkCollision(Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, Collision &colision){
    // distances to intersections
    std::vector<double> dist_intersections;

    walker.previous_location = walker.location;

    //cout << "-------------------" << endl; ;

    std::vector<Sphere*> spheres;

    find_all_intersections(walker, step, step_lenght + barrier_tickness, dist_intersections, spheres);
    
    if (dist_intersections.empty()){
        //cout << "dist_intersections.empty()" << endl;
        if (walker.location == Walker::intra){
            bool isinside = isPosInsideGlialCell(walker.pos_v, EPS_VAL);
            if (!isinside){
                colision.col_location = Collision::outside;
                walker.in_obj_index = -1;
                walker.in_obj_type = -1;
                walker.location = Walker::extra;
                // isinside = isPosInsideGlialCell_verbose(walker.pos_v, EPS_VAL);
                // isinside = isPosInsideGlialCell_verbose(walker.last_pos_v, EPS_VAL);
                // assert(0);
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


        for (auto i = 0; i < idx.size(); ++i) {
            int index = idx[i];
            Eigen::Vector3d pos = walker.pos_v + dist_intersections[index]*step;
            bool isnearEdge = true;
            if (walker.location == Walker::intra) 
                isnearEdge = !isPosInsideGlialCell(pos, -EPS_VAL);

            if (isnearEdge ){

                if (dist_intersections[index] <= step_lenght+barrier_tickness ){
                    
                                
                    //cout << "distance to sphere wall :" << sph.minDistance(pos)  << endl;
                    colision.type = Collision::hit;
                    colision.colision_point = pos;
                    colision.obstacle_ind = id;
                    colision.t = dist_intersections[index];


                    //Normal point
                    Eigen::Vector3d normal = (colision.colision_point - spheres[index]->P).normalized();
                    // Elastic bounce
                    Eigen::Vector3d temp_step = step;
                    Eigen::Vector3d ray =  (-colision.t*temp_step).normalized();//
                    double rn = ray.dot(normal);
                    temp_step = -ray + 2.0*normal*rn;

                    colision.bounced_direction = temp_step.normalized();
                    colision.perm_crossing = 0.;

                    
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
                        // cout << " outside coll 5 " << rn << " " << colision.t << " " << spheres[index]->id << " " << spheres[index]->object_id << endl;
                        // isPosInsideGlialCell_verbose(pos, EPS_VAL);    
                        // isPosInsideGlialCell_verbose(walker.last_pos_v, EPS_VAL);    
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

                    colision.perm_crossing = 0.;
                    return true;
                }
                else{
                    colision.type = Collision::null;
                    return false;
                }
            }
            
        }
        colision.type = Collision::null;

        return false;
        
    }
}

bool Glial::FindSphereinGlial(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<Sphere*> &sph){
  
    if (isNearGlialCell(position, distance_to_be_inside))
    {
        std::vector<std::vector<int>> spheres_id_to_check;
        
        int count = 0;
        // check soma
        for (auto axis = 0; axis < 3; ++axis) 
        {
            spheres_id_to_check.push_back(checkAxisForCollision_soma(position,distance_to_be_inside, axis));
            if(spheres_id_to_check[spheres_id_to_check.size()-1].size() > 0)
                ++count;
            
        }
        if(count == 3)
        {
            // find common ids in all 3 axes
            std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
            
            if (spheres_to_check_all_axes.size() > 0)
            {
                // soma
                if (soma.minDistance(position) <= distance_to_be_inside)
                    sph.push_back(&soma);
            }
        }
        
        spheres_id_to_check.clear();
        // check dendrites
        std::vector<int> dendrite_ids;
        bool isNearDendrite_ = isNearDendrite(position, distance_to_be_inside, dendrite_ids);
        if(isNearDendrite_)
        {
            for (auto d_id = 0; d_id < dendrite_ids.size(); ++d_id) 
            {
                count = 0;
                for (auto axis = 0; axis < 3; ++axis)
                {
                    spheres_id_to_check.push_back(checkAxisForCollision_dendrite(position,distance_to_be_inside, axis, dendrite_ids[d_id])); // check for collision along 1 axis
                    if(spheres_id_to_check[spheres_id_to_check.size()-1].size() > 0)
                        ++count;
                }

                 if(count == 3)
                {
                    // find common ids in all 3 axes
                    std::vector<int> spheres_to_check_all_axes = findCommonIntegers(spheres_id_to_check[0], spheres_id_to_check[1], spheres_id_to_check[2]);
                    for (auto i = 0; i < spheres_to_check_all_axes.size(); ++i) 
                    {
                        Sphere* sphere_to_check = &processes[dendrite_ids[d_id]][spheres_to_check_all_axes[i]];
                        if (sphere_to_check->minDistance(position) <= distance_to_be_inside)
                            sph.push_back(sphere_to_check);
                    }
                    spheres_to_check_all_axes.clear();
                }
                spheres_id_to_check.clear();
            }
        }
    }
    else
        return false;


    if (sph.size() > 0)
    {
        return true;
    }
    else
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
    for(unsigned d_id= 0 ; d_id < processes.size();d_id++)
    {   for(unsigned sph_id= 0 ; sph_id < processes[d_id].size();sph_id++)
        {
            processes[d_id][sph_id].percolation = this->percolation;

            if(processes[d_id][sph_id].percolation > 0.0){
                processes[d_id][sph_id].prob_cross_e_i = this->prob_cross_e_i;
                processes[d_id][sph_id].prob_cross_i_e = this->prob_cross_i_e;

                processes[d_id][sph_id].diffusivity_e = this->diffusivity_e;
                processes[d_id][sph_id].diffusivity_i = this->diffusivity_i;
                
            }
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
