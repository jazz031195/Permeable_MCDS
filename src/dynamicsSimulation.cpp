/** Main class of the dynamic Simulation.
 *  Everything related with the diffusion and collision of the particles is either implemented
 *  in this or used in this class.
 *
 * startSimulation method carries most of the computational burden.
 *
 * Project MC/DC Simulator
 * @author Jonathan
 * @version 1.42 stable
 */


#include "dynamicsSimulation.h"
#include "constants.h"
#include <algorithm>
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <chrono>     
#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <fstream>
#include "constants.h"
#include "collisionsphere.h"
#include "simerrno.h"
#include <chrono>
#include "simulablesequence.h"

using namespace Eigen;
using namespace std;
using namespace sentinels;

/**
 * DynamicsSimulation implementation
 */
DynamicsSimulation::DynamicsSimulation() {
    params.num_walkers = 1; //N
    params.num_steps   = 1; //T
    params.traj_file   = "";
    params.hit_file     = "";
    ini_pos_file       = "";
    completed          = 0;
    max_simulation_time = 0;

    tot_nbr_bounces = 0;
    tot_nbr_legal_crossings = 0;
    inside_bounces = 0;
    outside_bounces = 0;

    params.sim_duration = 1; //secs

    num_simulated_walkers = 0;

    params.diffusivity_intra = DIFF_CONST;
    params.diffusivity_extra = DIFF_CONST;
    
    step_lenght_intra = sqrt(6.0*params.diffusivity_intra*params.sim_duration/params.num_steps); // m
    step_lenght_extra = sqrt(6.0*params.diffusivity_extra*params.sim_duration/params.num_steps); // m

    //cout << "step_lenght_intra: " << step_lenght_intra << endl;
    //cout << "step_lenght_extra: " << step_lenght_extra << endl;
    //assert(0);

    params.write_traj = trajectory.write_traj = false;
    params.write_hit = trajectory.write_hit = false;
    params.write_txt = trajectory.write_txt   = false;

    nbr_walker_extra_ini = 0;
    nbr_walker_intra_ini = 0;
    nbr_walker_glials_ini = 0;
    nbr_walker_axons_ini = 0;
    nbr_walker_extra_final = 0;
    nbr_walker_intra_final = 0;
    nbr_walker_glials_final = 0;
    nbr_walker_axons_final = 0;
    nbr_walker_blood_vessels_ini = 0;
    nbr_walker_blood_vessels_final = 0;

    if(params.seed > 0){
        mt.seed(ulong(params.seed));
    }
    else{
        //Random seed
        std::random_device rd;
        mt.seed(rd());
    }

    print_expected_time = 1;
    icvf=0;
    intra_tries=0;
    total_tries=0;
    step_nbr = 0;
}

/**
 * @param configuration file
 */
DynamicsSimulation::DynamicsSimulation(std::string conf_file) {
    completed = 0;
    readConfigurationFile(conf_file);

    trajectory.initTrajectory(params);

    if(params.seed > 0){
        mt.seed(ulong(params.seed));
    }
    else{
        //Random seed
        std::random_device rd;
        mt.seed(rd());
    }
    print_expected_time = 1;

    this->step(1)=0;
    icvf=0;
    intra_tries=0;
    total_tries=0;
    step_nbr = 0;

    nbr_walker_extra_ini = 0;
    nbr_walker_intra_ini = 0;
    nbr_walker_glials_ini = 0;
    nbr_walker_axons_ini = 0;
    nbr_walker_extra_final = 0;
    nbr_walker_intra_final = 0;
    nbr_walker_glials_final = 0;
    nbr_walker_axons_final = 0;
    nbr_walker_blood_vessels_ini = 0;
    nbr_walker_blood_vessels_final = 0;

    tot_nbr_bounces = 0;
    tot_nbr_legal_crossings = 0;
    inside_bounces = 0;
    outside_bounces = 0;

}

/**
 * @param Parameter instance
 */
DynamicsSimulation::DynamicsSimulation(Parameters& params_) {
    params = params_;
    completed = 0;
    trajectory.initTrajectory(params);

    if(params.seed > 0){
        mt.seed(ulong(params.seed));
    }
    else{
        //Random seed
        std::random_device rd;
        mt.seed(rd());
    }

    print_expected_time  = 1;
    this->step(1)=0;
    icvf=0;
    intra_tries=0;
    total_tries=0;
    step_nbr = 0;

    nbr_walker_extra_ini = 0;
    nbr_walker_intra_ini = 0;
    nbr_walker_glials_ini = 0;
    nbr_walker_axons_ini = 0;
    nbr_walker_extra_final = 0;
    nbr_walker_intra_final = 0;
    nbr_walker_glials_final = 0;
    nbr_walker_axons_final = 0;
    nbr_walker_blood_vessels_ini = 0;
    nbr_walker_blood_vessels_final = 0;

    tot_nbr_bounces = 0;
    tot_nbr_legal_crossings = 0;
    inside_bounces = 0;
    outside_bounces = 0;
}

void DynamicsSimulation::initObstacleInformation(){


    if(params.collision_sphere_distance<= 0){
        params.collision_sphere_distance = inner_col_dist_factor;
    }
    double prob_cross_i_e, prob_cross_e_i, dse, dsi;

    //Cylinders list of index initialization

    for(unsigned i= 0 ; i < cylinders_list.size();i++){
        cylinders_deque.push_back(i);

        if(cylinders_list[i].percolation > 0.0){
            dse = sqrt(step_length_pref*cylinders_list[i].diffusivity_e);
            dsi = sqrt(step_length_pref*cylinders_list[i].diffusivity_i);

            prob_cross_i_e = cylinders_list[i].percolation * dsi * 2. / 3. / cylinders_list[i].diffusivity_i;
            prob_cross_e_i = cylinders_list[i].percolation * dse * 2. / 3. / cylinders_list[i].diffusivity_e; 

            cylinders_list[i].prob_cross_e_i = prob_cross_e_i / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
            cylinders_list[i].prob_cross_i_e = prob_cross_i_e / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
        }
    }

    walker.collision_sphere_cylinders.collision_list        = &cylinders_deque;
    walker.collision_sphere_cylinders.list_size             = unsigned(cylinders_deque.size());
    walker.collision_sphere_cylinders.big_sphere_list_end   = walker.collision_sphere_cylinders.list_size;
    
    //Axons list of index initialization
    int nbr_prints = 0;
    double null_prob = 0.0;
    for(unsigned i= 0 ; i < axons_list.size();i++){
        axons_deque.push_back(i);

        if(axons_list[i].percolation > 0.0){

            dse = sqrt(step_length_pref*axons_list[i].diffusivity_e);
            dsi = sqrt(step_length_pref*axons_list[i].diffusivity_i);

            prob_cross_i_e = axons_list[i].percolation * dsi * 2. / 3. / axons_list[i].diffusivity_i;
            prob_cross_e_i = axons_list[i].percolation * dse * 2. / 3. / axons_list[i].diffusivity_e; 
            axons_list[i].prob_cross_e_i = prob_cross_e_i / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
            axons_list[i].prob_cross_i_e = prob_cross_i_e / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
 
        }
    
    }

    for(unsigned i= 0 ; i < inner_axons_list.size();i++){

        if(inner_axons_list[i].percolation > 0.0){

            dse = sqrt(step_length_pref*inner_axons_list[i].diffusivity_e);
            dsi = sqrt(step_length_pref*inner_axons_list[i].diffusivity_i);

            prob_cross_i_e = inner_axons_list[i].percolation * dsi * 2. / 3. / inner_axons_list[i].diffusivity_i;
            prob_cross_e_i = inner_axons_list[i].percolation * dse * 2. / 3. / inner_axons_list[i].diffusivity_e; 
            inner_axons_list[i].prob_cross_e_i = prob_cross_e_i / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
            inner_axons_list[i].prob_cross_i_e = prob_cross_i_e / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e)); 
        }
    }

    walker.axons_collision_sphere.collision_list        = &axons_deque;
    walker.axons_collision_sphere.list_size             = unsigned(axons_deque.size());
    walker.axons_collision_sphere.big_sphere_list_end   = walker.axons_collision_sphere.list_size;

    //Glial cells list of index initialization

    for(unsigned i= 0 ; i < glials_list.size();i++){
        glials_deque.push_back(i);

        if(glials_list[i].percolation > 0.0){

            dse = sqrt(step_length_pref*glials_list[i].diffusivity_e);
            dsi = sqrt(step_length_pref*glials_list[i].diffusivity_i);
            
            prob_cross_i_e = glials_list[i].percolation * dsi * 2. / 3. / glials_list[i].diffusivity_i;
            prob_cross_e_i = glials_list[i].percolation * dse * 2. / 3. / glials_list[i].diffusivity_e; 
            glials_list[i].prob_cross_e_i = prob_cross_e_i / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
            glials_list[i].prob_cross_i_e = prob_cross_i_e / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
        }
    }

    walker.collision_sphere_glials.collision_list        = &glials_deque;
    walker.collision_sphere_glials.list_size             = unsigned(glials_deque.size());
    walker.collision_sphere_glials.big_sphere_list_end   = walker.collision_sphere_glials.list_size;


    // blood vessels list of index initialization
    for (unsigned i= 0 ; i < blood_vessels_list.size();i++){
        blood_vessels_deque.push_back(i);
    }
    walker.collision_sphere_blood_vessels.collision_list        = &blood_vessels_deque;
    walker.collision_sphere_blood_vessels.list_size             = unsigned(blood_vessels_deque.size());
    walker.collision_sphere_blood_vessels.big_sphere_list_end   = walker.collision_sphere_blood_vessels.list_size;


    // PLY index list initialization
    for(unsigned i= 0 ; i < plyObstacles_list.size();i++){
        std::vector<unsigned> jkr;
        for(unsigned t =0; t < plyObstacles_list[i].face_number;t++){
            jkr.push_back(t);
        }

        if(plyObstacles_list[i].percolation > 0.0){
            dse = sqrt(step_length_pref*plyObstacles_list[i].diffusivity_e);
            dsi = sqrt(step_length_pref*plyObstacles_list[i].diffusivity_i);

            prob_cross_i_e = plyObstacles_list[i].percolation * dsi * 2. / 3. / plyObstacles_list[i].diffusivity_i;
            prob_cross_e_i = plyObstacles_list[i].percolation * dse * 2. / 3. / plyObstacles_list[i].diffusivity_e; 

            plyObstacles_list[i].prob_cross_e_i = prob_cross_e_i / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
            plyObstacles_list[i].prob_cross_i_e = prob_cross_i_e / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));

            
        }

        ply_deque.push_back(jkr);
        walker.collision_sphere_ply.small_sphere_list_end.push_back(0);
        walker.collision_sphere_ply.big_sphere_list_end.push_back(unsigned(jkr.size()));
        walker.collision_sphere_ply.list_size++;
    }
    walker.collision_sphere_ply.collision_list = &ply_deque;

    //Spheres list of index initialization
    for(unsigned i= 0 ; i < spheres_list.size();i++){
        spheres_deque.push_back(i);   

        if(spheres_list[i].percolation > 0.0){

            double dse = sqrt(step_length_pref*spheres_list[i].diffusivity_e);
            double dsi = sqrt(step_length_pref*spheres_list[i].diffusivity_i);
            
            double prob_cross_i_e = spheres_list[i].percolation * dsi * 2. / 3. / spheres_list[i].diffusivity_i;
            double prob_cross_e_i = spheres_list[i].percolation * dse * 2. / 3. / spheres_list[i].diffusivity_e; 

            spheres_list[i].prob_cross_e_i = prob_cross_e_i / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
            spheres_list[i].prob_cross_i_e = prob_cross_i_e / (1.+ 0.5 * (prob_cross_e_i + prob_cross_i_e));
        }
    }
    walker.collision_sphere_spheres.collision_list        = &spheres_deque;
    walker.collision_sphere_spheres.list_size             = unsigned(spheres_deque.size());
    walker.collision_sphere_spheres.big_sphere_list_end   = walker.collision_sphere_spheres.list_size;
}

void DynamicsSimulation::updatePropagator(Eigen::Matrix3Xd& log_pos_r)
{
    for ( uint t =0; t <  params.record_prop_times.size(); t++){
        auto time = params.record_prop_times[t];

        Eigen::Vector3f displacement;

        displacement[0] = float(log_pos_r(0,0) - log_pos_r(0,time));
        displacement[1] = float(log_pos_r(1,0) - log_pos_r(1,time));
        displacement[2] = float(log_pos_r(2,0) - log_pos_r(2,time));

        for (uint i=0; i < params.prop_dirs.size();i++){

            Eigen::Vector3f direction = params.prop_dirs[i];

            Eigen::Vector3f d_projection = direction.dot(displacement)*direction;

            propagator.propagator_log[t][i] +=  d_projection.squaredNorm();
        }
    }
}

void DynamicsSimulation::normalizePropagator(float num_samples)
{

    if(num_samples<0)
        return;

    for (uint i =0 ; i < this->propagator.num_times; i++){
        for (uint j = 0 ; j < propagator.num_dirs; j++){
            propagator.propagator_log[i][j]/=num_samples;
        }
    }
}

void DynamicsSimulation::computeICVF()
{
    icvf = float(intra_tries)/float(total_tries);
}

bool DynamicsSimulation::finalPositionCheck()
{   
    int object_id, object_type;
    if(plyObstacles_list.size()>0 and sentinela.deport_illegals and params.obstacle_permeability <=0 and params.axon_obstacle_permeability <=0 and params.glial_obstacle_permeability <=0){
 
        bool isIntra = isInIntra(this->walker.pos_v, object_id, object_type);

        if((isIntra and this->walker.initial_location==Walker::extra) or ((!isIntra and this->walker.initial_location==Walker::intra))){
            //cout << "Im working" << endl;
            return true;
        }
    }
    return false;
}

void DynamicsSimulation::writePropagator(std::string path)
{
    if(params.write_bin){
        ofstream bout;
        bout.open(path.c_str(), std::ofstream::binary);

        if(!bout){
            std::cout << "Cannot open " << path << std::endl;
            return;
        }

        for (uint i =0 ; i < this->propagator.num_times; i++){
            for (uint j = 0 ; j < propagator.num_dirs; j++){
                float p = propagator.propagator_log[i][j];
                bout.write(reinterpret_cast<char *>(&p), sizeof(p));
            }
        }
        bout.close();
    }

    if(params.write_txt){

        ofstream tout;
        tout.open(path, std::ios::out);

        if(!tout){
            //TODO: Error handling
            std::cout << "Cannot open " << path.c_str()<< std::endl;
            return;
        }

        for (uint i =0 ; i < this->propagator.num_times; i++){
            for (uint j = 0 ; j < propagator.num_dirs; j++){
                float p = propagator.propagator_log[i][j];
                tout << p << " ";
            }
            tout << endl;
        }

        tout.close();
    }
}




void DynamicsSimulation::initSimulation()
{

    // Writes the header file and opens .traj file (if required)
    trajectory.initTrajWriter();
    trajectory.initHitWriter();
    trajectory.initFullWriter();
    
    // Initialize the walker trajectory log
    walker.setNumberOfSteps(params.num_steps);


    //time step = dt/T
    time_step = params.sim_duration/double(params.num_steps);
    time_dt =    0;
    last_time_dt=0;

    // Initial step length = sqrt(6*D*dt/T)
    step_length_pref  = 6.0*time_step;
    step_lenght_intra = sqrt(step_length_pref*params.diffusivity_intra);
    step_lenght_extra = sqrt(step_length_pref*params.diffusivity_extra);

    //cout << "step_lenght_intra: " << step_lenght_intra << endl;
    //cout << "step_lenght_extra: " << step_lenght_extra << endl;
    //assert(0);

    curr_step_lenght  = step_lenght_intra;
    curr_diffusivity  = params.diffusivity_intra;

    if (params.mix_compartment_diffusivities) {

        params.t_ex= params.t_ex*1e3;

        r1 = (1-params.f)/params.t_ex;
        r2 = params.f/params.t_ex;

        alpha = 1-exp(-r1*time_step);
        beta  = 1-exp(-r2*time_step);

    }

    //cout << " time_step = " << time_step << endl;

    //to predict time
    time(&start);

    //Initialize reading file for the walker position (if passed as parameter)
    std::ios_base::iostate exceptionMask = iniPos.exceptions() | std::ios::failbit;
    iniPos.exceptions(exceptionMask);

    // if a ini_walkers_file was passed
    if(params.ini_walkers_file.length() > 3){
        try {
            iniPos.open(params.ini_walkers_file);
            ini_pos_file_ini_index = (uint(id)*params.num_walkers)%params.ini_walkers_file_count;

            //We move the file start until we reach the initial position (MULTICORE SUPPPORT)
            for(unsigned i=0; i < ini_pos_file_ini_index; i++){
                double x,y,z;
                iniPos >> x; iniPos >> y; iniPos >> z;
            }
        }
        catch (std::ios_base::failure& e) {
            std::cerr << "Sim: " << id << " " << e.what() << '\n' << "Error loading walker initial positions. \n";
        }
    }

    initObstacleInformation(); 

    //Flags for the crossing and stuck particles. (numerical error sentinels)
    sentinela.deport_illegals = params.discard_illegals;
    sentinela.discard_stucks  = params.discard_stucks;

    if(params.log_propagator){
        propagator.num_dirs = uint(params.prop_dirs.size());
        propagator.num_times = uint(params.record_prop_times.size());
        propagator.initPropagator();
    }


    if(params.custom_sampling_area == false and voxels_list.size()>0){
        for(auto i = 0; i<3;i++){
           params.min_sampling_area[i]=voxels_list[0].min_limits[i];
           params.max_sampling_area[i]=voxels_list[0].max_limits[i];
        }
    }
}


bool DynamicsSimulation::expectedTimeAndMaxTimeCheck(unsigned w)
{
    /******  Some code to get a predicted time  **********/
    double completed_perc = int((w+1.0)/(params.num_walkers+1.0)*100.0);
    //unsigned tent_num_walkers = params.num_walkers/10 + 1;

    time(&now);

    second_passed = now - start;

    if(print_expected_time){
        if (completed_perc >= completed){

            if(this->completed <= 0.0){
                SimErrno::expectedTime(to_string(int(completed)), "Unknown",std::cout,true,"?","");
                std::cout.flush();
                completed+=10;
            }
            else if(completed_perc >= completed){
                std::cout << string(50,' ');
                std::cout << string(200,'\b');
                std::cout.flush();
                int steps_per_second = int(float(w*params.num_steps)/second_passed)*params.num_proc;

                if(steps_per_second > 0){
                string s_p =  std::to_string(steps_per_second);

                SimErrno::expectedTime(to_string(int(completed)), secondsToMinutes( (double(params.num_walkers - w+1))*second_passed/double(w+1)),
                                       std::cout,true,s_p,"");

                }
                else{
                    SimErrno::expectedTime(to_string(int(completed)), "Unknown",std::cout,true,"?","");
                    std::cout.flush();
                }
                completed = max(completed+10.0,completed_perc);
                std::cout.flush();
            }
        }
        else if( w == params.num_walkers-1){
            std::cout << string(50,' ');
            std::cout << string(200,'\b');
            std::cout.flush();
            SimErrno::expectedTime("100", "0 seconds",std::cout,true,"\n");
            std::cout.flush();
        }
    }
    /****** END Some code to get a predicted time  **********/


    if( (params.max_simulation_time>0) && (second_passed >= params.max_simulation_time) ){
        return true;
    }

    return false;
}

void DynamicsSimulation::writeDWSignal(SimulableSequence* dataSynth)
{
    //Writes the output data
    if(dataSynth){
        dataSynth->writeResultingData(params.output_base_name);
        if(params.log_phase_shift)
            dataSynth->writePhaseShiftDistribution(params.output_base_name);
    }
}

void DynamicsSimulation::iniWalkerPosition()
{
    //cout << "Initializing walker position..." << endl;
    walker.initial_location = Walker::unknown;
    walker.location         = Walker::unknown;
    walker.intra_extra_consensus = walker.intra_coll_count = walker.extra_coll_count=walker.rejection_count=0;
    int object_id = -1;
    int object_type = -1;
    //If the number of positions is less than the walkers, it restarts.
    if(iniPos.is_open()){
        double x,y,z;

        iniPos >> x; iniPos >> y; iniPos >> z;
        walker.setInitialPosition(x,y,z);

        if(++ini_pos_file_ini_index >= params.ini_walkers_file_count ){
            iniPos.clear();
            iniPos.seekg(0);
            ini_pos_file_ini_index = 0;
        }
    }
    else if (params.ini_delta_pos.size() > 0){
        walker.setRandomInitialPosition(Vector3d(double(params.ini_delta_pos[0]),double(params.ini_delta_pos[1]),double(params.ini_delta_pos[2])),
                Vector3d(double(params.ini_delta_pos[0]),double(params.ini_delta_pos[1]),double(params.ini_delta_pos[2])));
    }
    else if(params.ini_walker_flag.compare("intra")== 0){
        Vector3d intra_pos;
        getAnIntraCellularPosition(intra_pos, object_id, object_type);
        if (object_type == 0){
            nbr_walker_axons_ini++;
        }
        else if (object_type == 1){
            nbr_walker_glials_ini++;
        }
        walker.setInitialPosition(intra_pos);
        walker.intra_extra_consensus--;
        walker.initial_location = Walker::intra;
        walker.location = Walker::intra;
        walker.previous_location = Walker::intra;
        walker.in_obj_index = object_id;
        walker.in_obj_type = object_type;
        nbr_walker_intra_ini++;

    }
    else if(params.ini_walker_flag.compare("extra")== 0){
        Vector3d extra_pos;
        getAnExtraCellularPosition(extra_pos);
        walker.setInitialPosition(extra_pos);
        walker.intra_extra_consensus++;
        walker.initial_location = Walker::extra;
        walker.previous_location = Walker::extra;
        walker.location = Walker::extra;
        walker.in_obj_index = -1;
        walker.in_obj_type = -1;
        nbr_walker_extra_ini++;
    }
    else if(voxels_list.size() > 0 or params.custom_sampling_area){


        //cout << "params.max_sampling_area :" << params.max_sampling_area << endl;
        walker.setRandomInitialPosition(params.min_sampling_area,params.max_sampling_area);
        // Walker initial position - Required for multiple diffusivities
        bool isintra = isInIntra(walker.ini_pos, object_id, object_type, -barrier_tickness);
        bool isextra = isInExtra(walker.ini_pos, barrier_tickness);

        bool illegal = (!isintra && !isextra) || (isintra && isextra);

        while (illegal){
            walker.setRandomInitialPosition(params.min_sampling_area,params.max_sampling_area);
            isintra = isInIntra(walker.ini_pos, object_id, object_type, -barrier_tickness);
            isextra = isInExtra(walker.ini_pos, barrier_tickness);

            illegal = (!isintra && !isextra) || (isintra && isextra);
        }

        if (isintra){
            walker.initial_location = Walker::intra;
            walker.location = Walker::intra;
            walker.previous_location = Walker::intra;
            walker.in_obj_index = object_id;
            walker.in_obj_type = object_type;
            nbr_walker_intra_ini++;
            if (object_type == 0){
                nbr_walker_axons_ini++;
            }
            else if (object_type == 1){
                nbr_walker_glials_ini++;
            }
        }
        else {
            walker.initial_location = Walker::extra;
            walker.location = Walker::extra;
            walker.previous_location = Walker::extra;
            walker.in_obj_index = -1;
            walker.in_obj_type = -1;
            nbr_walker_extra_ini++;
        }

        if(params.computeVolume){
            isInIntra(walker.ini_pos, object_id, object_type, 0.0);
        }
        
    }
    else{
        walker.setInitialPosition(Vector3d(0,0,0));
    }
}


void DynamicsSimulation::initWalkerObstacleIndexes()
{
    // The outer collision sphere has a radius r = l*T 
    float outer_col_dist_factor = float(params.num_steps*curr_step_lenght);
    
    walker.initial_sphere_pos_v = walker.pos_v;
    walker.collision_sphere_cylinders.setBigSphereSize(outer_col_dist_factor);
    
    // The inner collision sphere has radius l*T*collision_sphere_distance
    //float inner_col_dist_factor = curr_step_lenght*sqrt(params.num_steps)*params.collision_sphere_distance;
    float inner_col_dist_factor = curr_step_lenght*5;
    
    walker.collision_sphere_cylinders.setSmallSphereSize(inner_col_dist_factor);

    // New version Cylinders obstacle selection
    walker.collision_sphere_cylinders.small_sphere_list_end = 0;
    walker.collision_sphere_cylinders.big_sphere_list_end = unsigned(cylinders_deque.size());

    // We add and remove the cylinder indexes that are or not inside sphere.
    for(unsigned i = 0 ; i < walker.collision_sphere_cylinders.list_size; i++ ){
        unsigned index = walker.collision_sphere_cylinders.collision_list->at(i);
        float dist = float(cylinders_list[index].minDistance(walker));
        if (dist < walker.collision_sphere_cylinders.small_sphere_distance){
            walker.collision_sphere_cylinders.pushToSmallSphere(i);
        }
    }
    
    // Axons Collision Sphere 
    walker.axons_collision_sphere.setBigSphereSize(outer_col_dist_factor);
    walker.axons_collision_sphere.setSmallSphereSize(inner_col_dist_factor);
    // New version Axons obstacle selection
    walker.axons_collision_sphere.small_sphere_list_end = 0;
    walker.axons_collision_sphere.big_sphere_list_end = unsigned(axons_deque.size());

    for (unsigned i = 0 ; i < walker.axons_collision_sphere.list_size; i++ ){
        unsigned index = walker.axons_collision_sphere.collision_list->at(i);
        float dist = float(axons_list[index].minDistance(walker));
        if (dist < walker.axons_collision_sphere.small_sphere_distance){
            walker.axons_collision_sphere.pushToSmallSphere(i);

        }
    }


    //* Glial Collision Sphere *//
    walker.collision_sphere_glials.setBigSphereSize(outer_col_dist_factor);
    walker.collision_sphere_glials.setSmallSphereSize(inner_col_dist_factor);
    // New version Glial obstacle selection
    walker.collision_sphere_glials.small_sphere_list_end = 0;
    walker.collision_sphere_glials.big_sphere_list_end = unsigned(glials_deque.size());
    // We add and remove the glial indexes that are or not inside sphere.
    for(unsigned i = 0 ; i < walker.collision_sphere_glials.list_size; i++ ){
        unsigned index = walker.collision_sphere_glials.collision_list->at(i);
        float dist = float(glials_list[index].minDistance(walker));
        if (dist < walker.collision_sphere_glials.small_sphere_distance){
            walker.collision_sphere_glials.pushToSmallSphere(i);
        }
    }

    //* Blood Vessels Collision Sphere *//
    walker.collision_sphere_blood_vessels.setBigSphereSize(outer_col_dist_factor);
    walker.collision_sphere_blood_vessels.setSmallSphereSize(inner_col_dist_factor);
    // New version Blood Vessels obstacle selection
    walker.collision_sphere_blood_vessels.small_sphere_list_end = 0;
    walker.collision_sphere_blood_vessels.big_sphere_list_end = unsigned(blood_vessels_deque.size());
    // We add and remove the blood vessels indexes that are or not inside sphere.
    for(unsigned i = 0 ; i < walker.collision_sphere_blood_vessels.list_size; i++ ){
        unsigned index = walker.collision_sphere_blood_vessels.collision_list->at(i);
        float dist = float(blood_vessels_list[index].minDistance(walker));
        if (dist < walker.collision_sphere_blood_vessels.small_sphere_distance){
            walker.collision_sphere_blood_vessels.pushToSmallSphere(i);
        }
    }

    //* PLY Collision Sphere *//
    
    walker.collision_sphere_ply.setBigSphereSize(outer_col_dist_factor);
    walker.collision_sphere_ply.setSmallSphereSize(inner_col_dist_factor);

    for(unsigned i = 0 ; i < walker.collision_sphere_ply.list_size; i++ )
    {
        walker.collision_sphere_ply.small_sphere_list_end[i] = 0;
        walker.collision_sphere_ply.big_sphere_list_end[i] = plyObstacles_list[i].face_number;
        for(unsigned t = 0 ; t < plyObstacles_list[i].face_number; t++){

            unsigned index = walker.collision_sphere_ply.collision_list->at(i)[t];
            float dist = float(plyObstacles_list[i].minDistance(walker,index));     

            if (dist > walker.collision_sphere_ply.big_sphere_distance)
            {
                walker.collision_sphere_ply.popFromBigSphere(i,t);
            }

            if (dist < walker.collision_sphere_ply.small_sphere_distance)
            {
                walker.collision_sphere_ply.pushToSmallSphere(i,t);
            }
        }
    }

    /* Sphere  */
        // The outer collision sphere has a radius r = l*T 
    walker.collision_sphere_spheres.setBigSphereSize(outer_col_dist_factor);
    
    // The inner collision sphere has radius l*T*collision_sphere_distance
    walker.collision_sphere_spheres.setSmallSphereSize(inner_col_dist_factor);

    // New version  obstacle selection
    walker.collision_sphere_spheres.small_sphere_list_end = 0;
    walker.collision_sphere_spheres.big_sphere_list_end = unsigned(spheres_deque.size());

    // We add and remove the sphere indexes that are or not inside sphere.
    for(unsigned i = 0 ; i < walker.collision_sphere_spheres.list_size; i++ ){
        unsigned index = walker.collision_sphere_spheres.collision_list->at(i);
        float dist = float(spheres_list[index].minDistance(walker));
        if (dist < walker.collision_sphere_spheres.small_sphere_distance){
            walker.collision_sphere_spheres.pushToSmallSphere(i);
        }
    }
}

void DynamicsSimulation::updateStepLength(double &l, const int &t){

    if (params.mix_compartment_diffusivities){
        if (t == 0){
            std::uniform_real_distribution<double> udist(0,1);
            double p = udist(mt);

            if (p < params.f){
                curr_step_lenght = step_lenght_intra;
                curr_diffusivity = params.diffusivity_intra;
            }
            else {
                curr_step_lenght = step_lenght_extra;
                curr_diffusivity = params.diffusivity_extra;
            }
            l = curr_step_lenght;
        }

        else{

            if (curr_step_lenght == step_lenght_intra){
                double prob_intra_to_extra = alpha;
                std::uniform_real_distribution<double> udist(0,1);
                double p = udist(mt);
                if (p < prob_intra_to_extra){
                    curr_step_lenght = step_lenght_extra;
                    curr_diffusivity = params.diffusivity_extra;
                    //cout <<" Switching to extra-cellular step length, time : " <<t << endl;
                }
                else{
                    curr_step_lenght = step_lenght_intra;
                    curr_diffusivity = params.diffusivity_intra;
                }
            }
            else if (curr_step_lenght == step_lenght_extra){
                double prob_extra_to_intra = beta;
                std::uniform_real_distribution<double> udist(0,1);
                double p = udist(mt);
                if (p < prob_extra_to_intra){
                    curr_step_lenght = step_lenght_intra;
                    curr_diffusivity = params.diffusivity_intra;
                    //cout <<" Switching to intra-cellular step length, time : " <<t << endl;
                }
                else{
                    curr_step_lenght = step_lenght_extra;
                    curr_diffusivity = params.diffusivity_extra;
                }
            }
            l = curr_step_lenght;
        }
    }
    
    else{

        if (walker.location == Walker::intra && walker.in_obj_type != blood_obstacle_type){
            l                = step_lenght_intra;
            curr_step_lenght = step_lenght_intra;
            curr_diffusivity = params.diffusivity_intra;
            //cout << "step length: " << l << " at time " << t <<  "walker.in_obj_type : " << walker.in_obj_type << "but not"<< blood_obstacle_type<<endl;
        }
        else if (walker.location == Walker::intra && walker.in_obj_type == blood_obstacle_type && walker.status != Walker::bouncing){

            Blood_Vessel bv = blood_vessels_list[walker.in_obj_index];
            double v;
            Eigen::Vector3d direction_flow;
            bv.WalkerVelocity(walker, v, direction_flow);

            l = float(v*params.sim_duration/float(params.num_steps));
            //cout << "Blood vessel flow step length: " << l << " at time " << t << " velocity " << v << endl;
            curr_step_lenght = l;
            curr_diffusivity = 0.0;
        }

        else if (walker.location == Walker::extra){
            l                = step_lenght_extra;
            curr_step_lenght = step_lenght_extra;
            curr_diffusivity = params.diffusivity_extra;
        } 
    }
}

void DynamicsSimulation::updateCollitionSphere(unsigned t)
{
    float inner_ball_size = walker.collision_sphere_ply.small_sphere_distance;
    float outer_ball_size = walker.collision_sphere_ply.big_sphere_distance;

    float sphere_sqrd_displacement = float((walker.initial_sphere_pos_v-walker.pos_v).norm());

    if(sphere_sqrd_displacement  + float(curr_step_lenght)   > outer_ball_size){
        initWalkerObstacleIndexes();
    }
    else if(sphere_sqrd_displacement + float(curr_step_lenght) > inner_ball_size    ){
        updateWalkerObstacleIndexes(t);
    }
}

void DynamicsSimulation::getAnIntraCellularPosition(Vector3d &intra_pos, int &object_id, int &object_type)
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);


    if(axons_list.size() <=0 and inner_axons_list.size() <=0 and cylinders_list.size() <=0 and plyObstacles_list.size() <= 0 and spheres_list.size() <= 0 and glials_list.size() <= 0 and blood_vessels_list.size() <=0){
        SimErrno::error("Cannot initialize intra-axonal walkers within the given substrate.",std::cout);
        SimErrno::error("There's no defined intra-axonal compartment (missing obstacles?)",std::cout);
        assert(0);
    }
    if(voxels_list.size()<=0){
        SimErrno::error("Cannot initialize intra-cellular walkers within the given substrate, no voxel.",std::cout);
        assert(0);
    }

    unsigned count = 0;
    while(true){

        if(count > 10000000){
            SimErrno::error("Cannot initialize intra-axonal walkers within the given substrate",std::cout);
            SimErrno::error("Max. number of tries to find an intra-celular compartment reached",std::cout);
            assert(0);
        }

        double x = double(udist(gen));
        double y = double(udist(gen));
        double z = double(udist(gen));

        x = x*(params.min_sampling_area[0]) + ( 1.0-x)*params.max_sampling_area[0];
        y = y*(params.min_sampling_area[1]) + ( 1.0-y)*params.max_sampling_area[1];
        z = z*(params.min_sampling_area[2]) + ( 1.0-z)*params.max_sampling_area[2];


       // std::cout << initialization_gap[2] << endl;
        Vector3d pos_temp = {x,y,z};

        std::vector<int> expected_object_types;
        if (inner_axons_list.size() > 0){
            expected_object_types.push_back(axon_obstacle_type);
        }
        else if (axons_list.size() > 0){
            expected_object_types.push_back(axon_obstacle_type);
        }
        if (glials_list.size() > 0){
            expected_object_types.push_back(glial_obstacle_type);
        }
        if (cylinders_list.size() > 0){
            expected_object_types.push_back(cyl_obstacle_type);
        }
        if (blood_vessels_list.size() > 0){
            expected_object_types.push_back(blood_obstacle_type);
        }

        bool isintra = isInIntra(pos_temp, object_id, object_type, -barrier_tickness);

        if(checkIfPosInsideVoxel(pos_temp) && isintra && (std::find(expected_object_types.begin(), expected_object_types.end(), object_type) != expected_object_types.end()) ){
            
            intra_pos = pos_temp;
            walker.initial_location = Walker::intra;
            walker.location = Walker::intra;
            walker.previous_location = Walker::intra;
            walker.in_obj_index = object_id;
            walker.in_obj_type = object_type;
            
            return;
        }
        count++;
    }
}

void DynamicsSimulation::getAnExtraCellularPosition(Vector3d &extra_pos)
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);



    if(voxels_list.size()<=0){
        SimErrno::error("Cannot initialize extra-cellular walkers within the given substrate, no voxel.",std::cout);
        assert(0);
    }

    unsigned count = 0;

    while(true){

        if(count > 1000000){
            SimErrno::error("Cannot initialize extra-cellular walkers within the given substrate",std::cout);
            SimErrno::error("Max. number of tries to find an extra-celular compartment reached",std::cout);
            assert(0);
        }

        double x = double(udist(gen));
        double y = double(udist(gen));
        double z = double(udist(gen));

        x = x*(params.min_sampling_area[0]) + ( 1.0-x)*params.max_sampling_area[0];
        y = y*(params.min_sampling_area[1]) + ( 1.0-y)*params.max_sampling_area[1];
        z = z*(params.min_sampling_area[2]) + ( 1.0-z)*params.max_sampling_area[2];

        Vector3d pos_temp = {x,y,z};

        bool isextra = isInExtra(pos_temp, barrier_tickness);

        if(checkIfPosInsideVoxel(pos_temp) && isextra){
            extra_pos = pos_temp;
            walker.initial_location = Walker::extra;
            walker.location = Walker::extra;
            walker.previous_location = Walker::extra;
            walker.in_obj_index = -1;
            walker.in_obj_type = -1;
            return;
        }
        count++;
    }
}

bool DynamicsSimulation::checkIfPosInsideVoxel(Vector3d &pos)
{

    for(unsigned v = 0; v < voxels_list.size();v++){
        if(     pos[0] - voxels_list[v].min_limits[0] > barrier_tickness &&
                pos[1] - voxels_list[v].min_limits[1] > barrier_tickness &&
                pos[2] - voxels_list[v].min_limits[2] > barrier_tickness &&
                pos[0] - voxels_list[v].max_limits[0] < -barrier_tickness &&
                pos[1] - voxels_list[v].max_limits[1] < -barrier_tickness &&
                pos[2] - voxels_list[v].max_limits[2] < -barrier_tickness)
            return true;
    }

    return false;
}

//TODO: Use t to decrease the size of the sphere.
void DynamicsSimulation::updateWalkerObstacleIndexes(unsigned t_)
{

    //float outher_col_dist_factor = float(params.num_steps-t_+1.0*curr_step_lenght);
    float outher_col_dist_factor = float(params.num_steps-t_+1.0*step_lenght_intra);
    
    walker.collision_sphere_ply.setBigSphereSize(outher_col_dist_factor);

    walker.initial_sphere_pos_v = walker.pos_v;

    //Cylinders obstacle update.
    walker.collision_sphere_cylinders.small_sphere_list_end = 0;

    for(unsigned i = 0 ; i < walker.collision_sphere_cylinders.big_sphere_list_end; i++ )
    {
        unsigned index = walker.collision_sphere_cylinders.collision_list->at(i);
        float dist    = float(cylinders_list[index].minDistance(walker));

        if (dist > walker.collision_sphere_cylinders.big_sphere_distance)
        {
            walker.collision_sphere_cylinders.popFromBigSphere(i);
        }
        if (dist < walker.collision_sphere_cylinders.small_sphere_distance)
        {
            walker.collision_sphere_cylinders.pushToSmallSphere(i);
        }
    }
    
    //Axons obstacle update.
    //cout << "updateWalkerObstacleIndexes" << endl;

    //auto time = std::chrono::high_resolution_clock::now();
    walker.axons_collision_sphere.small_sphere_list_end = 0;

    for(unsigned i = 0 ; i < walker.axons_collision_sphere.big_sphere_list_end; i++ )
    {
        unsigned index = walker.axons_collision_sphere.collision_list->at(i);
        double dist = axons_list[index].minDistance(walker);

        if (dist > walker.axons_collision_sphere.big_sphere_distance)
        {
            walker.axons_collision_sphere.popFromBigSphere(i);

        }
        if (dist < walker.axons_collision_sphere.small_sphere_distance)
        {
            walker.axons_collision_sphere.pushToSmallSphere(i);
        }
        
    }

    //auto time2 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( time2 - time ).count();
    //cout << "Time to update axons: " << duration << endl;
    //assert(0);
    //cout << "end updateWalkerObstacleIndexes" << endl;

    //PLY update obstacle
    for(unsigned i = 0 ; i < walker.collision_sphere_ply.list_size; i++ )
    {
        walker.collision_sphere_ply.small_sphere_list_end[i] = 0;

        for(unsigned t = 0 ; t < walker.collision_sphere_ply.big_sphere_list_end[i]; t++){

            unsigned triangle_index = walker.collision_sphere_ply.collision_list->at(i)[t];
            float dist = float(plyObstacles_list[i].minDistance(walker,triangle_index));

            if (dist > walker.collision_sphere_ply.big_sphere_distance)
            {
                walker.collision_sphere_ply.popFromBigSphere(i,t);
            }

            if (dist < walker.collision_sphere_ply.small_sphere_distance)
            {
                walker.collision_sphere_ply.pushToSmallSphere(i,t);
            }
        }
    }

    //Spheres obstacle update.
    walker.collision_sphere_spheres.small_sphere_list_end = 0;

    for(unsigned i = 0 ; i < walker.collision_sphere_spheres.big_sphere_list_end; i++ )
    {
        unsigned index = walker.collision_sphere_spheres.collision_list->at(i);
        float dist    = float(spheres_list[index].minDistance(walker));

        if (dist > walker.collision_sphere_spheres.big_sphere_distance)
        {
            walker.collision_sphere_spheres.popFromBigSphere(i);
        }
        if (dist < walker.collision_sphere_spheres.small_sphere_distance)
        {
            walker.collision_sphere_spheres.pushToSmallSphere(i);
        }
    }

}

string DynamicsSimulation::secondsToMinutes(double t)
{
    if(t < 60){
        return std::to_string(int(t)) + " seconds";
    }

    int mins    = int(t/60.0);

    int seconds = int(t - mins*60);

    string min_s =  (mins>1)?" minutes":" minute";
    return std::to_string(mins) +  min_s +  " and " + std::to_string(seconds) + " seconds";

}


bool DynamicsSimulation::isInsideCylinders(Vector3d &position, int &object_id, double distance_to_be_inside)
{
    Walker tmp;
    tmp.setInitialPosition(position);

    //track the number of positions checks for intra/extra positions

    for(unsigned i = 0 ; i < cylinders_list.size(); i++){

        double dis = cylinders_list[i].minDistance(tmp);

        if( dis <= distance_to_be_inside ){
            intra_tries++;
            object_id = i;
            return true;
        }
    }
    object_id = -1;
    return false;
}

bool DynamicsSimulation::isInsideAxons(Eigen::Vector3d &position, int &object_id, double distance_to_be_inside)
{
    if (inner_axons_list.size() > 0) {
        for (unsigned i = 0; i < inner_axons_list.size() ; i++){
            bool isinside = inner_axons_list[i].isPosInsideObstacle(position,  distance_to_be_inside);
            if (isinside){
                object_id = i;
                return true;
            }
        }
    }
    else{
        for (unsigned i = 0; i < axons_list.size() ; i++){
    
            bool isinside = axons_list[i].isPosInsideObstacle(position,  distance_to_be_inside);
            if (isinside){
                object_id = i;
                return true;
            }
        }
    }
    object_id = -1;
    return false;
}

bool DynamicsSimulation::isOutsideAxons(Eigen::Vector3d &position, int &object_id, double distance_to_be_inside)
{

    for (unsigned i = 0; i < axons_list.size() ; i++){
 
        bool isinside = axons_list[i].isPosInsideObstacle(position,  distance_to_be_inside);
        if (isinside){
            object_id = i;
            return false;
        }
    }
    object_id = -1;
    return true;
}

bool DynamicsSimulation::isOutsideCylinders(Eigen::Vector3d &position, int &object_id, double distance_to_be_inside)
{
    for (unsigned i = 0; i < cylinders_list.size() ; i++){
 
        bool isinside = cylinders_list[i].minDistance(position) < distance_to_be_inside;
        if (isinside){
            object_id = i;
            return false;
        }
    }
    object_id = -1;
    return true;
}


bool DynamicsSimulation::isInsidePLY(Eigen::Vector3d &position, double distance_to_be_inside)
{

    //1) We find the closest PLY and triangle bases on the triangle's center
    Walker tmp;
    tmp.setInitialPosition(position);

    double t,min_t = 1e6;
    unsigned min_j_index = 0;
    int min_i_index = -1;
    for (unsigned i=0; i < plyObstacles_list.size(); i++){
        for (unsigned j=0; j < plyObstacles_list[i].face_number; j++){
            t = (position - plyObstacles_list[i].faces[j].center).squaredNorm();
            if(t< min_t){
                min_i_index = i;
                min_j_index = j;
                min_t = t;
            }
        }
    }

    //2) We corroborate by casting an infinite ray and checking collisions

    Eigen::Vector3d ray = (-position + plyObstacles_list[min_i_index].faces[min_j_index].center).normalized();
    Collision collision_temp;

    double new_min_t = 1e6;
    for (unsigned i=0; i < plyObstacles_list.size(); i++){
        for (unsigned j=0; j < plyObstacles_list[i].face_number; j++){
            plyObstacles_list[i].faces[j].stepIntersects_MT(tmp,ray,1e8,collision_temp);

            if(collision_temp.type == Collision::hit and new_min_t > collision_temp.t){
                new_min_t = collision_temp.t;
                min_i_index = i;
                min_j_index = j;
            }
        }
    }

    //3) Finally we check the sign of the closest collision. The sign indicates either intra or extra.
    if(min_i_index >= 0){
        Eigen::Vector3d normal;
        plyObstacles_list[min_i_index].faces[min_j_index].getNormal(normal);

        //Orientation respect the triangle
        double dot = ((position - plyObstacles_list[min_i_index].faces[min_j_index].center).normalized()).dot(normal);
        if (dot < distance_to_be_inside){
            intra_tries++;
            return true;
        }
    }

    return false;
}


bool DynamicsSimulation::isInsideSpheres(Vector3d &position, double distance_to_be_inside)
{
    Walker tmp;
    tmp.setInitialPosition(position);

    //track the number of positions checks for intra/extra positions

    for(unsigned i = 0 ; i < spheres_list.size(); i++){

        double dis = spheres_list[i].minDistance(tmp);

        if( dis <= distance_to_be_inside ){
            intra_tries++;
            return true;
        }
    }

    return false;
}

bool DynamicsSimulation::isInsideGlials(Eigen::Vector3d &position, int &object_id, const double& distance_to_be_inside)
{

    for (unsigned i = 0; i < glials_list.size() ; i++){
        bool isinside = glials_list[i].isPosInsideObstacle(position,  distance_to_be_inside);
        if (isinside){
            object_id = int(i);

            return true;
        }
    }
    object_id = -1;

    return false;
}

bool DynamicsSimulation::isInsideBloodVessels(Eigen::Vector3d &position, int &object_id, const double& distance_to_be_inside)
{
    
    for (unsigned i = 0; i < blood_vessels_list.size() ; i++){

        bool isinside = blood_vessels_list[i].isPosInsideObstacle(position,  distance_to_be_inside);
        if (isinside){
            object_id = int(i);
            return true;
        }
    }
    object_id = -1;

    return false;
}

bool DynamicsSimulation::isInIntra(Vector3d &position, int &object_id, int& object_type, double distance_to_be_intra_ply)
{

    bool isIntra = false;
    total_tries++;
    object_type = -1;
    int ax_id, glial_id, cyl_id, bv_id;
    bool isinside_axons = false, isinside_glial = false, isinside_cyl = false, isinside_bv = false;

    if(cylinders_list.size()>0){
        isinside_cyl= this->isInsideCylinders(position, cyl_id, distance_to_be_intra_ply);
        isIntra|= isinside_cyl;
        if (isinside_cyl){
            object_type = cyl_obstacle_type;
        }

    }

    if (inner_axons_list.size()>0){
        isinside_axons = this->isInsideAxons(position, ax_id, distance_to_be_intra_ply);
        isIntra|= isinside_axons;
        if (isinside_axons){
            object_type = axon_obstacle_type;
        }
    }

    if(glials_list.size()>0){
        isinside_glial = this->isInsideGlials(position, glial_id, distance_to_be_intra_ply);
        isIntra|= isinside_glial;
        if (isinside_glial){
            object_type = glial_obstacle_type;

        }
    }

    if (blood_vessels_list.size()>0){
        isinside_bv = this->isInsideBloodVessels(position, bv_id, distance_to_be_intra_ply);
        isIntra|= isinside_bv;
        if (isinside_bv){
            object_type = blood_obstacle_type;

        }
    }

    if(plyObstacles_list.size()>0){
        isIntra|=isInsidePLY(position,distance_to_be_intra_ply);
        assert(0);
    }

    if(spheres_list.size()>0){
        isIntra|= this->isInsideSpheres(position,distance_to_be_intra_ply);   

    }

    if (isinside_axons){
        object_id = ax_id;
    }
    else if (isinside_glial){
        object_id = glial_id;
    }
    else if (isinside_cyl){
        object_id = cyl_id;
    }
    else if (isinside_bv){
        object_id = bv_id;
    }
    else{
        object_id = -1;
    }
 

    return isIntra;
}


bool DynamicsSimulation::isInExtra(Eigen::Vector3d &position,  double distance_to_be_intra_ply)
{

    bool ok = true; int dummy;

    // cylinders/axons/spheres: ideally use the same distance-threshold idea for them too
    if (!cylinders_list.empty()) ok = ok && isOutsideCylinders(position, dummy, distance_to_be_intra_ply);
    if (!axons_list.empty())     ok = ok && isOutsideAxons(position, dummy, distance_to_be_intra_ply);
    if (!plyObstacles_list.empty()) ok = ok && !isInsidePLY(position, distance_to_be_intra_ply);
    if (!spheres_list.empty())      ok = ok && !isInsideSpheres(position, distance_to_be_intra_ply);
    if (!glials_list.empty())   ok = ok && !this->isInsideGlials(position, dummy, distance_to_be_intra_ply);
    if (!blood_vessels_list.empty())   ok = ok && !this->isInsideBloodVessels(position, dummy, distance_to_be_intra_ply);
    return ok;
}


void DynamicsSimulation::startSimulation(SimulableSequence *dataSynth) {

    //Initialize values, arrays and files.
    initSimulation();

    //Alias of the step length, may vary when the time step is dynamic. Default in extra-cellular space
    double l = curr_step_lenght;
    bool back_tracking;
    int num_crossed = 0;

    list_walkers_intra = std::vector<int>(params.num_steps,0);
    list_walkers_extra = std::vector<int>(params.num_steps,0);


    /*********************   WARNING  **********************/
    /*                                                     */
    /*                 DYNAMIC SIMULATION CORE             */
    /*                                                     */
    /*********************   WARNING  **********************/
    unsigned w=0;

    //cout << "params.num_walkers :" << params.num_walkers << endl;
    
    for (w = 0 ; w < params.num_walkers; w++)
    {

        // start timer 
        auto start = std::chrono::high_resolution_clock::now();
        //flag in case there was any error with the particle.
    
        back_tracking = false;

        cout << "Progress :" << w << "/" << params.num_walkers << "( " << double(w*100/params.num_walkers) << " %)" << endl;
        

        walker.setIndex(w);

        walker.normal = {0,0,0};
  
        // Initialize the walker initial position
        iniWalkerPosition();

        // Update step length based on the walker initial position in space
        updateStepLength(l, 0);

        // Selects only obstacles that are close enough to collide and the ones inside a collision sphere
        initWalkerObstacleIndexes();
        
        //Initial position;
        walker.setRealPosLog(walker.pos_r,0);
        walker.setVoxPosLog (walker.pos_v,0);

        for(unsigned t = 1 ; t <= params.num_steps; t++) //T+1 steps in total (avoid errors)
        {         

            //Get the time step in milliseconds         
            getTimeDt(last_time_dt,time_dt,l,dataSynth,t,time_step);

            //Generates a random oriented step of size l
            generateStep(step,l);

            // Moves the particle. Checks collision and handles bouncing.
            
            try{

                //updateWalkerPosition(step);
                updateWalkerPosition(step, t);

            }
            catch(Sentinel::ErrorCases error){

                // Possible errors, or numerical un-handed cases should end here.
                sentinela.deportationProcess(walker,w,t,back_tracking,params,id);

                if ( (error == Sentinel::ErrorCases::stuck) || (error == Sentinel::ErrorCases::crossed)){
                    break;
                }

                if ( error == Sentinel::rejected  )
                    continue;
            }     


            // Saves the final particle position after bouncing in the time t.
            walker.setRealPosLog(walker.pos_r,t);
            walker.setVoxPosLog (walker.pos_v,t);

            // Save the collision 
            walker.setColision(walker.collision_in, walker.collision_ext, walker.crossing_in, walker.crossing_ext, t);           

            // Update step length and current diffusivity based on the walker position in space
            updateStepLength(l, t);

            if (walker.location != Walker::intra){
                updateCollitionSphere(t);
            }

            walker.steps_count++;
            walker.rejection_count = 0;
            
            if (walker.location == Walker::intra){
                list_walkers_intra[t-1]+= 1;
            }
            else{
                list_walkers_extra[t-1]+= 1;
            }

        }// end for t

        // 2. Stop the timer
        auto stop = std::chrono::high_resolution_clock::now();

        // 3. Calculate the duration
        // You can change std::chrono::milliseconds to std::chrono::seconds or std::chrono::microseconds
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        // 4. Print the result
        //std::cout << "Time taken for all steps: " << duration.count() << " milliseconds" << std::endl;

        start = std::chrono::high_resolution_clock::now();

        
        //if(!back_tracking)
        //    if(finalPositionCheck()){
        //        back_tracking=true;
        //        w--;
        //    }
        
        //If there was an error, we don't compute the signal or write anything.
        if(back_tracking){
            continue;
        }

        if (walker.location == Walker::intra){
            nbr_walker_intra_final ++;
            if (walker.in_obj_type == axon_obstacle_type){
                nbr_walker_axons_final ++;
            }
            else if (walker.in_obj_type == glial_obstacle_type){
                nbr_walker_glials_final ++;
            }
        }
        else{
            nbr_walker_extra_final ++;
        }

        dataSynth->update_phase_shift(this->time_step,walker.pos_r_log);
        dataSynth->update_DWI_signal(walker);
        //Write the positions.
        if (params.write_location || params.write_txt){
            trajectory.writePosition(walker, w);
        }
        else{
            trajectory.writePosition(walker.pos_r_log, walker.collision_in_log, walker.collision_ext_log, walker.crossing_in_log, walker.crossing_ext_log);
        }
        
        if(params.log_propagator){
            //Update Propagator
            updatePropagator(walker.pos_r_log);
        }

        //Displays the remained expected time and check for the time limit.
        if(expectedTimeAndMaxTimeCheck(w)){
            std::cout << "\n" << SH_BG_LIGHT_YELLOW <<  "[Warning]" << SH_DEFAULT << "  Sim: " << id << " "
                 << "Max time limit reached: Simulation halted after "<< ++w << " spins" << endl;
            break;
        }

        // 2. Stop the timer
        stop = std::chrono::high_resolution_clock::now();

        // 3. Calculate the duration
        // You can change std::chrono::milliseconds to std::chrono::seconds or std::chrono::microseconds
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        // 4. Print the result
        //std::cout << "Time taken fo signal: " << duration.count() << " milliseconds" << std::endl;




    }// for w
   

    
    

    /*********************   WARNING  **********************/
    /*                                                     */
    /*         END OF THE DYNAMIC SIMULATION CORE          */
    /*                                                     */
    /*********************   WARNING  **********************/

    num_simulated_walkers = w;

    if(num_simulated_walkers<= params.num_walkers){

        trajectory.reWriteHeaderFile(num_simulated_walkers);

        this->params.num_walkers = num_simulated_walkers;
    }

    if(params.log_propagator){
        normalizePropagator(num_simulated_walkers);
    }

    //computes the ICVF from the initialization.
    computeICVF();

    // Info display.
    time(&now);
    second_passed = difftime(now,start);
    if(params.verbatim)
        SimErrno::info(" Sim: " + to_string(id) + " Simulation ended after: "
                       + secondsToMinutes(second_passed) + " seconds" ,std::cout);


    // Writes the final DWI signal, and the phase shift.
    if(params.log_opp)
        writeDWSignal(dataSynth);

    return;
}

DynamicsSimulation::~DynamicsSimulation() {
    if(iniPos.is_open())
        iniPos.close();
}

/**
 * @param conf_file_path
 * @return void
 */
void DynamicsSimulation::readConfigurationFile(std::string conf_file_path) {
    params.readSchemeFile(conf_file_path);
}

/**
 * @return void
 */
void DynamicsSimulation::generateStep(Vector3d & step, double l) {

    if(walker.status == Walker::on_object){
        step = walker.next_direction.normalized();
        return;
    }

    if(walker.in_obj_type ==blood_obstacle_type && walker.location == Walker::intra && walker.status != Walker::on_voxel && walker.status != Walker::bouncing){
        Blood_Vessel bv = blood_vessels_list[walker.in_obj_index];
        double v;
        Eigen::Vector3d direction_flow;
        bv.WalkerVelocity(walker, v, direction_flow);
        step = direction_flow;
        step.normalize();
        return;
    }

    std::uniform_real_distribution<double> dist(0,1);

    /* Unbiased random direction*/
    double theta  = 2.0*M_PI*dist(mt);
    double cosPhi = 2.0*dist(mt)-1.0;
    double cosTh  = cos(theta);
    double sinTh  = sin(theta);
    double sinPhi = sqrt(1.0-cosPhi*cosPhi);

    step(0) = l*cosTh*sinPhi;
    step(1) = l*sinTh*sinPhi;
    step(2) = l*cosPhi;

    step.normalize();

}

void DynamicsSimulation::generateDirectedStep(Vector3d &new_step, Vector3d &direction){

    std::uniform_real_distribution<double> dist(0,1);

    /* Unbiased random direction*/
    double theta  = 2.0*M_PI*dist(mt);
    double cosPhi = 2.0*dist(mt)-1.0;
    double cosTh  = cos(theta);
    double sinTh  = sin(theta);
    double sinPhi = sqrt(1.0-cosPhi*cosPhi);

    new_step(0) = cosTh*sinPhi;
    new_step(1) = sinTh*sinPhi;
    new_step(2) = cosPhi;

    new_step.normalize();

    double rn = direction.dot(new_step);

    if(rn<0.0)
        new_step*=-1.0;
}

/**
 * @return True if a bouncing in needed.
 */
//bool DynamicsSimulation::updateWalkerPosition(Eigen::Vector3d& step) {
bool DynamicsSimulation::updateWalkerPosition(Eigen::Vector3d& step, unsigned &t) {

    //cout <<"updateWalkerPosition" << endl;
  
    //new step to take
    Vector3d bounced_step = step.normalized(),end_point;
    Vector3d previous_real_position, previous_voxel_position, real_pos, voxel_pos;

    //Save the walker initial position.
    walker.getVoxelPosition(previous_voxel_position);
    walker.getRealPosition(previous_real_position);

    // Collision instance to save manage the collision (in Spanish).
    Collision collision;

    // True when the particle needs to be bounced and updates.
    bool bounced = false;
    bool update_walker_status  = false;

    // Maximum displacement. Is updated after each bouncing (if any)
    double tmax =  curr_step_lenght;

    // Clears the status of the sentinel.
    sentinela.clear();

    // Reset collision state of walker
    walker.collision_in = 0;
    walker.collision_ext = 0;
    walker.crossing_in = 0;
    walker.crossing_ext = 0;

    unsigned bouncing_count = 0;



    do{
        //cout <<"---------------------" << endl;
        bounced = false;
        bouncing_count++;
        walker.is_allowed_to_cross = false;
        walker.previous_location = walker.location;
        
        // Checks the number of bouncing per step.
        walker.steps_count++;
        // True if there was a collision and the particle needs to be bounced.
        update_walker_status |= checkObstacleCollision(bounced_step, tmax, end_point, collision);

        if (std::isnan(collision.t)) {
            cout << "collision.t: " << collision.t << endl;
            assert(0);
        }

        // Updates the position and bouncing direction.
        if(update_walker_status){


            //bounced = updateWalkerPositionAndHandleBouncing(bounced_step,tmax,collision);
            bounced = updateWalkerPositionAndHandleBouncing(bounced_step,tmax,collision, t);

            Eigen::Vector3d p0  = walker.pos_v;
            if (std::isnan(p0.x()) || std::isnan(p0.y()) || std::isnan(p0.z())) {
                cout << "\n[FATAL ERROR] Walker position is NaN!" << endl;
                cout << "collision.type: " << collision.type << endl;
                cout << "collision.t: " << collision.t << endl;
                assert(0);
            }

            // restarts the variables.
            
            update_walker_status = false;
            collision.type = Collision::null;
            collision.col_location  = Collision::unknown;
            collision.perm_crossing = 0.0;
            collision.t = INFINITY_VALUE;
            
        }
        else{
            if (collision.type == Collision::null){
                // clear from previous status
                walker.status = Walker::free;
                walker.next_direction = {0,0,0};
            }
        }

        //cout << "walker.previous_location: " << walker.previous_location << endl;
        //cout << "walker.location: " << walker.location << endl;
        //cout << "walker.object_id: " << walker.in_obj_index << endl;

        sentinela.checkErrors(walker,params,(plyObstacles_list.size() == 0),bouncing_count);

    }while(bounced);

    if(tmax >= 0.0){

        // Update the walker position after the bouncing (or not)
        walker.getRealPosition(real_pos);
        Eigen::Vector3d adapted_step;

        if (walker.normal == Eigen::Vector3d {0,0,0}){
            walker.setRealPosition(real_pos  + tmax*bounced_step);
        } 
        else{
            adapted_step =  findMirrorStep(bounced_step, walker.normal);
            walker.setRealPosition(real_pos  + tmax*adapted_step);
        } 


        walker.getVoxelPosition(voxel_pos);
        walker.setVoxelPosition(voxel_pos+ tmax*bounced_step);

    }


    return false;
}

bool DynamicsSimulation::checkObstacleCollision(Vector3d &bounced_step,double &tmax, Eigen::Vector3d& end_point,Collision& collision)
{

    Collision collision_tmp;
    collision_tmp.type = Collision::null;
    collision_tmp.t = INFINITY_VALUE;

    //Origin O
    Eigen::Vector3d ray_origin;
    walker.getVoxelPosition(ray_origin);

    if (std::isnan(ray_origin.x()) || std::isnan(ray_origin.y()) || std::isnan(ray_origin.z())) {
        cout << "\n[FATAL ERROR] ray_origin is NaN!" << endl;
        assert(0);
    }

    //To keep track of the closest collision
    double max_collision_distance = tmax;


    // The collision checks the five  possible obstacles in this order: Voxel, Cylinders, Axons, PLY, spheres.
    // The closest collision is kept at the end.
    
    //Check Voxel limits
    for(unsigned int i = 0 ; i < voxels_list.size(); i++ )
    {
        //voxels_list[i].CheckCollision(walker,bounced_step,tmax,collision_tmp);
        voxels_list[i].CheckCollision(walker,bounced_step,tmax,collision_tmp);
        handleCollisions(collision,collision_tmp,max_collision_distance,i);
    }

    if ((cylinders_list).size()>0){
        for(unsigned int i = 0 ; i < walker.collision_sphere_cylinders.small_sphere_list_end; i++ )
        {
            unsigned index = walker.collision_sphere_cylinders.collision_list->at(i);
            cylinders_list[index].checkCollision(walker,bounced_step,tmax,collision_tmp);
            handleCollisions(collision,collision_tmp,max_collision_distance,index);
        }
    }

    //For each Axon Obstacle
    if ((axons_list).size()>0 || (inner_axons_list).size()>0){
        // intra walkers


        if (walker.location== Walker::intra ){

            if (walker.in_obj_type == axon_obstacle_type && walker.in_obj_index != -1){
                
                if (inner_axons_list.size() > 0) {
    
                    (inner_axons_list)[walker.in_obj_index].checkCollision(walker,bounced_step,tmax,collision_tmp);
                }
                else{
                    (axons_list)[walker.in_obj_index].checkCollision(walker,bounced_step,tmax,collision_tmp);
                }
                handleCollisions(collision,collision_tmp,max_collision_distance,walker.in_obj_index);   

            }
        }
        // extra walkers or unknown
        else {
            for(unsigned int i = 0 ; i < walker.axons_collision_sphere.small_sphere_list_end; i++ ){
            //for (unsigned int i = 0 ; i < (axons_list).size(); i++ ){
                unsigned index = walker.axons_collision_sphere.collision_list->at(i);
                //unsigned index = i;
                (axons_list)[index].checkCollision(walker,bounced_step,tmax,collision_tmp);
                handleCollisions(collision,collision_tmp,max_collision_distance,index);  
            }
        }
    }
    
    //For each Glial Obstacle
    if ((glials_list).size()>0 ){
        
        // intra walkers
        if (walker.location== Walker::intra ){
            
            if (walker.in_obj_type == glial_obstacle_type && walker.in_obj_index != -1){
                (glials_list)[walker.in_obj_index].checkCollision(walker,bounced_step,tmax,collision_tmp);
                handleCollisions(collision,collision_tmp,max_collision_distance,walker.in_obj_index);  
            } 
        }
        // extra walkers or unknown
        else {

            for(unsigned int i = 0 ; i < walker.collision_sphere_glials.small_sphere_list_end; i++ ){
            //for (unsigned int i = 0 ; i < (glials_list).size(); i++ ){
                unsigned index = walker.collision_sphere_glials.collision_list->at(i);
                //unsigned index = i;
                (glials_list)[index].checkCollision(walker,bounced_step,tmax,collision_tmp);
                handleCollisions(collision,collision_tmp,max_collision_distance,index);     
            }
        }
        
    }

    // for each blood vessel obstacle
    if ((blood_vessels_list).size()>0 ){
        
        // intra walkers
        if (walker.location== Walker::intra ){
            
            if (walker.in_obj_type == blood_obstacle_type && walker.in_obj_index != -1){
               (blood_vessels_list)[walker.in_obj_index].checkCollision(walker,bounced_step,tmax,collision_tmp);
                handleCollisions(collision,collision_tmp,max_collision_distance,walker.in_obj_index);  
            } 
        }
        // extra walkers or unknown
        else {

            for(unsigned int i = 0 ; i < walker.collision_sphere_blood_vessels.small_sphere_list_end; i++ ){
            //for (unsigned int i = 0 ; i < (blood_vessels_list).size(); i++ ){
                unsigned index = walker.collision_sphere_blood_vessels.collision_list->at(i);
                //unsigned index = i;
                (blood_vessels_list)[index].checkCollision(walker,bounced_step,tmax,collision_tmp);
                handleCollisions(collision,collision_tmp,max_collision_distance,index);     
            }
        }
        
    }


    //For each PLY Obstacles
    for(unsigned int i = 0 ; i < walker.collision_sphere_ply.collision_list->size(); i++ )
    {

        plyObstacles_list[i].checkCollision(walker,bounced_step,tmax,collision_tmp, walker.collision_sphere_ply.collision_list->at(i),
                                            walker.collision_sphere_ply.small_sphere_list_end[i]);

        handleCollisions(collision,collision_tmp,max_collision_distance,i);
    }

    //For each Sphere Obstacles
    for(unsigned int i = 0 ; i < walker.collision_sphere_spheres.small_sphere_list_end; i++ )
    {
        unsigned index = walker.collision_sphere_spheres.collision_list->at(i);
        spheres_list[index].checkCollision(walker,bounced_step,tmax,collision_tmp);
        handleCollisions(collision,collision_tmp,max_collision_distance,index);
    }

    step_nbr+= 1;

    return collision.type != Collision::null;
}


void DynamicsSimulation::handleCollisions(Collision &collision, Collision &collision_2, double &max_collision_distance, unsigned indx)
{
    // nothing to do;
    if (collision_2.type == Collision::null)
        return;

    collision_2.obstacle_ind = int(indx);

    if (collision.type == Collision::hit || collision.type == Collision::boundary){
        if(collision_2.doIHaveMorePiorityThan(collision)){
            collision = collision_2;
            max_collision_distance = collision_2.t;
            collision.obstacle_ind = int(indx);
        }
        return;
    }

    if(collision.type == Collision::near ){
        if (collision_2.type == Collision::hit || collision_2.type == Collision::boundary){
            collision = collision_2;
            max_collision_distance = collision_2.t;
            collision.obstacle_ind = int(indx);
        }
        return;
    }

    // if we get here means that collision.type = 'null'
    if(collision_2.type == Collision::near){

        collision = collision_2;
        collision.obstacle_ind = int(indx);

        return;
    }

    collision = collision_2;
}


void DynamicsSimulation::mapWalkerIntoVoxel(const Eigen::Vector3d &bounced_step, const Collision &collision,const double &barrier_thickness)
{
  
    walker.setRealPosition(walker.pos_r + collision.t*bounced_step);

    Eigen::Vector3d voxel_pos = walker.pos_v + (collision.t)*bounced_step;

    bool mapped = false;
    for(int i = 0 ; i < 3; i++)
    {
        if ( fabs(voxel_pos[i] -  voxels_list[0].min_limits[i]) <= EPS_VAL){
            voxel_pos[i] = voxels_list[0].max_limits[i];
            mapped = true;
        }
        else if ( fabs(voxel_pos[i] - voxels_list[0].max_limits[i]) <= EPS_VAL){
            voxel_pos[i] = voxels_list[0].min_limits[i];
            mapped = true;
        }
    }

    walker.setVoxelPosition(voxel_pos);


    if (mapped){
        initWalkerObstacleIndexes();
    }
}
/*
void DynamicsSimulation::mapWalkerIntoVoxel_tortuous(Eigen::Vector3d& bounced_step, Collision &collision)
{
    //cout << "MAPPING " << endl;
    walker.setRealPosition(walker.pos_r + collision.t*bounced_step);
    Eigen::Vector3d position;
    //cout << "mapped " << endl;
    if (walker.location== Walker::extra){
        getAnExtraCellularPosition(position);
    }
    else{
        int ax_id;
        getAnIntraCellularPosition(position, ax_id);
        walker.in_ax_index = ax_id;
    }
    walker.setVoxelPosition(position);
    initWalkerObstacleIndexes();
    

}
*/
bool DynamicsSimulation::elasticBounceAgainstVoxel(const Eigen::Vector3d& previous_pos, const Eigen::Vector3d& current_pos, Eigen::Vector3d &normal, Eigen::Vector3d &step){
    
    Eigen::Vector3d ray = (-step).normalized();
    bool hit_wall = false;
    normal = {0, 0, 0};

    // Check each axis independently WITH A DIRECTIONAL GUARD
    for(int i = 0 ; i < 3; i++) {
        
        // Hit MIN limit AND moving outwards (negative direction)?
        if (current_pos[i] <= voxels_list[0].min_limits[i] + EPS_VAL && step[i] < 0.0) {
            normal[i] = 1.0; 
            hit_wall = true;
        }
        // Hit MAX limit AND moving outwards (positive direction)?
        else if (current_pos[i] >= voxels_list[0].max_limits[i] - EPS_VAL && step[i] > 0.0) {
            normal[i] = 1.0; 
            hit_wall = true;
        }
    }

    if (hit_wall) {
        normal.normalize();
        double rn = ray.dot(normal);
        step = -ray + 2.0 * normal * rn;
        normal = {abs(normal[0]), abs(normal[1]), abs(normal[2])};
        return true;
    } 
    
    return false; 
}

Eigen::Vector3d DynamicsSimulation::findMirrorStep(const Eigen::Vector3d& bounced_step, const Eigen::Vector3d& normal){

    Eigen::Vector3d mirroredVector = bounced_step;
    if (normal == Eigen::Vector3d {0,0,0}){
        return bounced_step;
    }
    else{
        for (int i = 0; i < 3; ++i) {
            if (normal[i] != 0){
                mirroredVector[i] = -bounced_step[i];
            }
        }
    }

    return mirroredVector;
}


void DynamicsSimulation::mapWalkerIntoVoxel_tortuous(const Eigen::Vector3d& bounced_step, Collision &collision)
{
    Eigen::Vector3d previous_v_pos = walker.pos_v;

    walker.setVoxelPosition(walker.pos_v + collision.t*bounced_step); 
    if (walker.normal == Eigen::Vector3d {0,0,0}){
        walker.setRealPosition(walker.pos_r + collision.t*bounced_step);
    } 
    else{
        Eigen::Vector3d adapted_step =  findMirrorStep(bounced_step, walker.normal);
        walker.setRealPosition(walker.pos_r + collision.t*adapted_step);
    } 
    
    Eigen::Vector3d temp_step = bounced_step;
    Eigen::Vector3d normal = {0,0,0} ;
    if (walker.location == Walker::intra && walker.in_obj_type == blood_obstacle_type){

        bool hit_wall = false;
        
        // 1. Detect which wall we hit (with directional safety guards)
        for(int i = 0 ; i < 3; i++) {
            if (walker.pos_v[i] <= voxels_list[0].min_limits[i] + EPS_VAL && bounced_step[i] < 0.0) {
                normal[i] = 1.0; 
                hit_wall = true;
            } else if (walker.pos_v[i] >= voxels_list[0].max_limits[i] - EPS_VAL && bounced_step[i] > 0.0) {
                normal[i] = 1.0; 
                hit_wall = true;
            }
        }

        if (hit_wall) {

            // 2. The 180-degree reversal: Go straight backwards down the pipe
            collision.bounced_direction = -bounced_step;
            
            // 3. Toggle the mirror state for accurate pos_r tracking
            for(int i = 0; i < 3; ++i) {
                if (normal[i] != 0) {
                    walker.normal[i] = (walker.normal[i] == 0) ? 1 : 0;
                }
            }
 
            initWalkerObstacleIndexes();
        } else {
            // SAFETY CATCH: Ignore fake microscopic hits
            collision.bounced_direction = bounced_step;
        }
    }
    else{
        bool mapped = elasticBounceAgainstVoxel(previous_v_pos, walker.pos_v,normal, temp_step);

        if (mapped) { 
            collision.bounced_direction = temp_step.normalized();
            
            for(int i = 0; i < 3; ++i) {
                if (normal[i] != 0) {
                    walker.normal[i] = (walker.normal[i] == 0) ? 1 : 0;
                }
            }
            initWalkerObstacleIndexes();
        } else {
            // SAFETY CATCH: If it was a fake microscopic hit, just keep moving forward 
            // without toggling the mirror states!
            collision.bounced_direction = bounced_step;
        }

    }

    

} 

void DynamicsSimulation::getTimeDt(double &last_time_dt, double &time_dt, double &l, SimulableSequence* dataSynth, unsigned t, double time_step)
{
    last_time_dt = time_step*(t-1);
    time_dt = time_step*(t);
    
    if(dataSynth){
        if(dataSynth->dynamic){
            if (walker.in_obj_type ==blood_obstacle_type && walker.location == Walker::intra && walker.status != Walker::bouncing){
                Blood_Vessel bv = blood_vessels_list[walker.in_obj_index];
                double v;
                Eigen::Vector3d direction_flow;
                bv.WalkerVelocity(walker, v, direction_flow);
                last_time_dt = dataSynth->time_steps[t-1];
                time_dt = dataSynth->time_steps[t];
                l = v * (time_dt - last_time_dt);
            }
            else{
                last_time_dt = dataSynth->time_steps[t-1];
                time_dt = dataSynth->time_steps[t];
                l = sqrt(6.0*(curr_diffusivity*(time_dt - last_time_dt)));
            }
        }
    }
}

//bool DynamicsSimulation::updateWalkerPositionAndHandleBouncing(Vector3d &bounced_step, double &tmax, Collision &collision)
bool DynamicsSimulation::updateWalkerPositionAndHandleBouncing(Vector3d &bounced_step, double &tmax, Collision &collision, unsigned &t)
{
    if (std::isnan(collision.t)) {
        cout << "\n[FATAL ERROR] collision.t is NaN!" << endl;
        assert(0);
    }
    // To avoid numerical errors.
    double min_step_length =  barrier_tickness;

    Eigen::Vector3d real_pos, voxel_pos;
    walker.getRealPosition(real_pos);
    walker.getVoxelPosition(voxel_pos);

    bool bounced = false;

    //Sets the status of the walker
    walker.status = Walker::free;

    if (tmax <= min_step_length)
    {

        tmax = 0.0;
        return false;   
    }    

    if (collision.type == Collision::leaked) {
        if (walker.location == Walker::intra){
            walker.intra_extra_consensus++;
            walker.location = Walker::extra;
            walker.in_obj_index = -1;
            walker.in_obj_type = -1;
        }
        else if (walker.location == Walker::extra){
            walker.intra_extra_consensus--;
            walker.location = Walker::intra;
            walker.in_obj_index = -1;
            walker.in_obj_type = -1;
        }
        return false;
    }

    if(collision.type == Collision::hit && collision.col_location != Collision::voxel)
    {
   
           //If the collision was really close we can't trust the normal direction;
        if(collision.t < 1e-10 && walker.status != walker.bouncing){
            sentinela.rejected_step = true;
            Eigen::Vector3d direction = -step;
            generateDirectedStep(walker.next_direction,direction);
            walker.status = Walker::on_object;
            return false;
        }
        tot_nbr_bounces++;

        // if we go from inside to outside
        if (collision.col_location == Collision::inside){
            inside_bounces++;
        }
        else if (collision.col_location == Collision::outside){
            outside_bounces++;
        }

        bounced = true;
        walker.status = Walker::bouncing;

        double displ = collision.t;

        // Membrane crossing due to permeability --> Tmax needs to be updated because diffusivity may have changed!
        int crossed = 0, col_loc = collision.col_location;

        if (collision.perm_crossing > EPS_VAL){
            
            crossed = 1;
            walker.is_allowed_to_cross = true;
            //cout << "walker.is_allowed_to_cross = true" << endl;
            tot_nbr_legal_crossings++;
            // crossing the membrane
            if(collision.col_location == Collision::inside){
                
                tmax = sqrt(params.diffusivity_extra/params.diffusivity_intra) * (tmax-displ);
                                
                walker.intra_extra_consensus--;
                walker.location = Walker::extra;
                //cout << "walker.location = Walker::extra (permeability)" << endl;
                walker.in_obj_type = -1;
                walker.in_obj_index = -1;
                /*
               std::cout << "Walker was in :" << walker.previous_location << " but is now in :" << walker.location << endl;
               std::cout << " walker is now extra after crossing" << endl;
               std::cout << "previous_location:" << walker.previous_location  << endl;
               std::cout << "location:" << walker.location  << endl;
                */
            // Save that walker hit the membrane - For validation purpose
                walker.crossing_in++;
                walker.collision_in++;

            }
            else if(collision.col_location == Collision::outside){
                
                tmax = sqrt(params.diffusivity_intra/params.diffusivity_extra) *(tmax -displ);
  
                walker.intra_extra_consensus++;
                walker.location = Walker::intra;
                //cout << "walker.location = Walker::intra (permeability)" << endl;
                walker.in_obj_index = collision.obstacle_ind;
                walker.in_obj_type = collision.obstacle_type;

                /*
               std::cout << "Walker was in :" << walker.previous_location << " but is now in :" << walker.location << endl;
               std::cout << " walker is now intra after crossing" << endl;
               std::cout << "previous_location:" << walker.previous_location  << endl;
               std::cout << "location:" << walker.location  << endl;
                */
                
            // Save that walker hit the membrane - For validation purpose
                walker.crossing_ext++;
                walker.collision_ext++;

            }
            else{
                tmax-=displ;
                //cout << "previous_location:" << walker.previous_location  << endl;
                //cout << "location:" << walker.location  << endl;

            }

            if(walker.initial_location == Walker::unknown){
                walker.initial_location = walker.location;
            }        

        }
        else{

            tmax -= displ;

            if(collision.col_location == Collision::inside){
                walker.location = Walker::intra;
                //cout << "walker.location = Walker::intra (no permeability)" << endl;
                walker.in_obj_index = collision.obstacle_ind;
                walker.in_obj_type = collision.obstacle_type;
                /*
               std::cout << "Walker was in :" << walker.previous_location << " but is now in :" << walker.location << endl;
               std::cout << " walker stays intra" << endl;
               std::cout << "previous_location:" << walker.previous_location  << endl;
               std::cout << "location:" << walker.location  << endl;
                */
            // Save that walker hit the membrane - For validation purpose
                walker.collision_in++;
            }
            if(collision.col_location == Collision::outside){
                walker.location = Walker::extra;
                //cout << "walker.location = Walker::extra (no permeability)" << endl;
                walker.in_obj_index = -1;
                walker.in_obj_type = -1;
                /*
               std::cout << "Walker was in :" << walker.previous_location << " but is now in :" << walker.location << endl;
               std::cout << " walker stays extra" << endl;
               std::cout << "previous_location:" << walker.previous_location  << endl;
               std::cout << "location:" << walker.location  << endl;
                */
            // Save that walker hit the membrane - For validation purpose
                walker.collision_ext++;
            }
            
        }
        
        if (tmax < 0.0){
            tmax = 0.0;
            bounced = false;
        }

        
        //We update the position.
        Eigen::Vector3d adapted_step;
        if (walker.normal == Eigen::Vector3d {0,0,0}){
            walker.setRealPosition(real_pos + displ*bounced_step);
        } 
        else{

            adapted_step =  findMirrorStep(bounced_step, walker.normal);
            walker.setRealPosition(real_pos + displ*adapted_step);
            
        } 
        walker.setVoxelPosition(voxel_pos  + displ*bounced_step);
        

        bounced_step = collision.bounced_direction;

        // Save collision - for validation purpose
        //trajectory.writeFullCollision(collision.collision_point, crossed, col_loc, t, walker.index);
    
    }
    else if(collision.type == Collision::hit && collision.col_location == Collision::voxel)
    {

        bounced = true;

        walker.status = Walker::on_voxel;

        mapWalkerIntoVoxel_tortuous(bounced_step,collision);

        tmax-=collision.t;

    }
    else if(collision.type == Collision::near){
        //sentinela.rejected_step   = true;
        Eigen::Vector3d direction = -bounced_step; //WARNING: deberiamos usar el bounced step.
        generateDirectedStep(walker.next_direction,direction);
        walker.status = Walker::on_object;

        if(collision.col_location == Collision::inside){
            walker.location = Walker::intra;
            walker.in_obj_index = collision.obstacle_ind;
            walker.in_obj_type = collision.obstacle_type;

        // Save that walker hit the membrane - For validation purpose
            walker.collision_in++;
        }
        else if(collision.col_location == Collision::outside){
            walker.location = Walker::extra;
            walker.in_obj_index = -1;
            walker.in_obj_type = -1;
        // Save that walker hit the membrane - For validation purpose
            walker.collision_ext++;
        }

        //Save info on collision - Validation purpose
        int crossed=0, col_loc = collision.col_location;
        //trajectory.writeFullCollision(collision.collision_point, crossed, col_loc, t, walker.index);
        
        return false;
    }
    else if(collision.type == Collision::degenerate){
        sentinela.rejected_step = true;
        Eigen::Vector3d direction = -step;
        generateDirectedStep(walker.next_direction,direction);
        walker.status = Walker::on_object;
        return false;
    }
   
    return bounced;
}



void DynamicsSimulation::setDuration(const double &duration)
{
    params.sim_duration = duration;
    trajectory.dyn_duration = duration;
}

void DynamicsSimulation::setWalkersNum(const unsigned &N)
{
    params.num_walkers  = N;
    trajectory.N = N;
}

void DynamicsSimulation::setStepsNum(const unsigned &T)
{
    params.num_steps = T;
    trajectory.T = T;
}


