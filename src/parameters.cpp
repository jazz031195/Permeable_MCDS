#include "parameters.h"
#include <fstream>
#include <iostream>
#include "constants.h"
#include "simerrno.h"
using namespace std;

Parameters::Parameters()
{
    //Dummy initial values;
    scale_from_stu      = false;
    seed                = -1;
    save_phase_shift    = false;
    write_traj          = false;
    write_hit           = false;
    write_full_c        = false;
    write_txt           = false;
    write_bin           =  true;

    hex_packing = false;
    hex_packing_radius      = 0;
    hex_packing_separation  = 0;

    gamma_packing   = false;
    gaussian_packing= false;
    neuron_packing  = false;
    uniform_packing = false;                      /*!< flag, true if a gamma distribution of spheres will be initialized        */
    hex_packing     = false;
    packing_cyl     = false;
    packing_s       = false;

    packing_icvf    = 0;
    packing_output_configuration = 0.0;

    // For compilation
    gamma_output_conf=0;
    gamma_icvf=EPS_VAL;
    gamma_output_configuration=0.0;
    gamma_num_cylinders=0;
    concentration = 0;

    gamma_num_axons=0;

    packing_output_conf = false;

    ini_walkers_file = "";
    num_proc    = 0;
    verbatim    = false;

    record_phase_times.clear();
    record_pos_times.clear();
    subdivisions_file = "";

    computeVolume = false;
    custom_sampling_area = false;
    for (auto i= 0;i<3; i++)
        min_sampling_area[i] = max_sampling_area[i] = 0.0;
}

void Parameters::readSchemeFile(std::string conf_file_path)
{

    ifstream in(conf_file_path);

    if(!in){
        cout << "[ERROR] Configuration file" << endl;
        return;
    }

    string tmp="";
    int number_lines = 0;
    while((in >> tmp) && (str_dist(tmp,"<END>") >= 2) ){
        //cout << "read line number :" << number_lines << endl;
        //cout << tmp << endl;

        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"n") == 0){
            in >> num_walkers;
        }
        else if(str_dist(tmp,"c") == 0){
            in >> concentration;
        }
        else if(str_dist(tmp,"t") == 0){
            in >> num_steps;
        }
        else if(str_dist(tmp,"duration") <= 1){
            in >> sim_duration;
        }
        else if(str_dist(tmp,"scheme_file") <= 1){
            in >> scheme_file;
        }
        else if(str_dist(tmp,"diffusivity_intra") <= 1){
            in >> diffusivity_intra;
        }
        else if(str_dist(tmp,"diffusivity_extra") <= 1){
            in >> diffusivity_extra;
        }
        else if( (str_dist(tmp,"out_traj_file_index") <= 2) or (str_dist(tmp,"exp_prefix") <= 2)) {
            in >> traj_file;
            hit_file= traj_file;
            output_base_name = traj_file;
        }
        else if( str_dist(tmp,"ini_walkers_file") <= 3){
            in >> ini_walkers_file;
        }
        else if(str_dist(tmp,"write_txt") <= 1){
            in >> write_txt;
        }
        else if(str_dist(tmp,"write_bin") <= 1){
            in >> write_bin;
        }
        else if(str_dist(tmp,"write_traj_file") <= 2){
            in >> write_traj;
        }
        else if(str_dist(tmp,"write_hit_file") <= 2){
            in >> write_hit;
        }
        else if(str_dist(tmp,"write_full_c_file") <= 2){
            in >> write_full_c;
        }
        else if(str_dist(tmp,"scale_from_stu") <= 2){
            in >> scale_from_stu;
        }
        else if(str_dist(tmp,"seed") <= 1){
            in >> seed;
        }
        else if(str_dist(tmp,"<obstacle>") == 0){
            readObstacles(in);
        }
        else if ((str_dist(tmp,"<voxels>") == 0) or (str_dist(tmp,"<voxel>") == 0)){
            readVoxels(in);
        }
        else if(str_dist(tmp,"num_process") <= 1 || str_dist(tmp,"processors") <= 1){
            in >> num_proc;
        }
        else if(str_dist(tmp,"<log>") == 0){
            readInfoGatheringParams(in);
        }
        else if(str_dist(tmp,"<sampling_area>") ==0 || str_dist(tmp,"<spawning_area>") ==0){
            float tmp;
            in >> tmp; min_sampling_area[0] = tmp;
            in >> tmp; min_sampling_area[1] = tmp;
            in >> tmp; min_sampling_area[2] = tmp;
            in >> tmp; max_sampling_area[0] = tmp;
            in >> tmp; max_sampling_area[1] = tmp;
            in >> tmp; max_sampling_area[2] = tmp;

            custom_sampling_area = true;
        }
        else if(str_dist(tmp,"<delta>") == 0)
        {
            float tmp;
            //Save the three positions
            in >> tmp;ini_delta_pos.push_back(tmp);
            in >> tmp;ini_delta_pos.push_back(tmp);
            in >> tmp;ini_delta_pos.push_back(tmp);
        }
        else if(str_dist(tmp,"ini_walkers_pos") <= 2)
        {
            in >> ini_walker_flag;
            std::transform(ini_walker_flag.begin(), ini_walker_flag.end(), ini_walker_flag.begin(), ::tolower);

            if(str_dist(ini_walker_flag,"intra") > 1 && str_dist(ini_walker_flag,"extra") > 1 && str_dist(ini_walker_flag,"delta") > 1){
                SimErrno::error("Wrong walker's initialization position (Did you mean: 'walker_ini_file'?) ",cout);

                assert(0);
            }
        }
        else if(str_dist(tmp,"subdivisions_number") <= 1)
        {
            in >> number_subdivisions;
            subdivision_flag |= (number_subdivisions>1)?true:false;
            if(number_subdivisions > 100){
                SimErrno::error("Unrealistic number of resulting subdivision voxels : " + std::to_string(number_subdivisions) + "^3",cout);

                assert(0);
            }
        }
        else if(str_dist(tmp,"subdivisions_file") <= 1)
        {
            in >> subdivisions_file;
            readSubdivisionFile();
            subdivision_flag |= (subdivisions.size()>0);
        }
        /*
        else if((str_dist(tmp,"permeability") <= 1) || (str_dist(tmp,"obstacle_permeability") <= 2))
        
        {
            readPermeability(in);
        }
        else if((str_dist(tmp,"<permeability>") <= 1))
        {
            readPermeability(in);
        }

        */
        else if( str_dist(tmp,"coll_sphere_dist") <= 2 || str_dist(tmp,"colision_dist") <= 2 || str_dist(tmp,"sphere_size") <= 2)
        {
            in >> collision_sphere_distance;
        }
        else if( str_dist(tmp,"verbatim") <= 2 )
        {
           this->verbatim = true;
        }
        else if( str_dist(tmp,"log_opp") <= 1 )
        {
           this->log_opp = true;
        }
        else if( str_dist(tmp,"compute_volume") <= 1 )
        {
           this->computeVolume = true;
        }
        else if( str_dist(tmp,"log_phase_shift") <= 2 )
        {
           this->log_phase_shift = true;
        }
        else if(str_dist(tmp,"max_sim_time") <= 2 || str_dist(tmp,"max_simulation_time") <= 3){
            in >> max_simulation_time;

            if(max_simulation_time<=0){
                SimErrno::error("Max simulation time must be greater than 1 second. ",cout);
                assert(0);
            }
        }
        else if( str_dist(tmp,"deportation") <= 2 )
        {
            in >> discard_illegals;
        }
        else if( str_dist(tmp,"discard_stucks") <= 2 )
        {
            in >> discard_stucks;
        }
        else{
            if( str_dist(tmp.substr(0,2),"</") > 0 )
                SimErrno::warning("Parameter: " + tmp + " Unknown",cout);
        }
        number_lines += 1;
    }

    if(scale_from_stu){
        //m^2/s to mm^2/ms
        diffusivity_intra*=m2_to_mm2/s_to_ms;
        diffusivity_extra*=m2_to_mm2/s_to_ms;
        //seconds to ms
        sim_duration*=s_to_ms;
    }
    cout << "number_subdivisions :" << number_subdivisions << endl;
    if(number_subdivisions>1){
        addSubdivisions();
    }

    in.close();
    return;
}

//Set Methods

void Parameters::setNumWalkers(unsigned N)
{
    num_walkers = N;
}

void Parameters::setNumSteps(unsigned T)
{
    num_steps = T;
}

void Parameters::setDiffusivity(double Di, double De)
{
    diffusivity_intra = Di;
    diffusivity_extra = De;
}

void Parameters::setSimDuration(double duration)
{
    sim_duration = duration;
}

void Parameters::setWriteTrajFlag(bool write_bin)
{
    write_traj = write_bin;
}
void Parameters::setWriteHitFlag(bool write_bin)
{
    write_hit = write_bin;
}

void Parameters::setWriteFullFlag(bool write_full_c_)
{
    write_full_c = write_full_c_;
}


void Parameters::setWriteTextFlag(bool write_txt_)
{
    write_txt = write_txt_;
}

void Parameters::setMinLimits(Eigen::Vector3d min_limits_)
{
    min_limits = min_limits_;

}

void Parameters::setMaxLimits(Eigen::Vector3d max_limits_)
{
    max_limits = max_limits_;
}

void Parameters::setTrajFileName(std::string traj_file_)
{
    traj_file = traj_file_;
}

void Parameters::setHitFileName(std::string hit_file_)
{
    hit_file = hit_file_;
}

void Parameters::setOutputBaseFileName(std::string output_base_name_)
{
    output_base_name = output_base_name_;
}

void Parameters::iniWalkersFileName(std::string ini_walkers_file_)
{
    ini_walkers_file = ini_walkers_file_;
}

void Parameters::setSchemeFileName(std::string scheme_file_)
{
    scheme_file = scheme_file_;
}


//GET METHODS

unsigned Parameters::getNumWalkers()
{
    return num_walkers;
}

unsigned Parameters::getNumSteps()
{
    return num_steps;
}

double Parameters::getDiffusivity_intra()
{
    return diffusivity_intra;
}

double Parameters::getDiffusivity_extra()
{
    return diffusivity_extra;
}
bool Parameters::getWriteTrajFlag()
{
    return write_traj;
}
bool Parameters::getWriteHitFlag()
{
    return write_hit;
}

bool Parameters::getWriteFullFlag()
{
    return write_full_c;
}

bool Parameters::getWriteTextFlag()
{
    return write_txt;
}

Eigen::Vector3d Parameters::getMinLimits()
{
    return min_limits;
}

Eigen::Vector3d Parameters::getMaxLimits()
{
    return max_limits;
}

std::string Parameters::getTrajFileName()
{
    return traj_file;
}

std::string Parameters::getHitFileName()
{
    return hit_file;
}

std::string Parameters::getOutputBaseFileName()
{
    return output_base_name;
}

std::string Parameters::getIniWalkersFileName()
{
    return ini_walkers_file;
}

std::string Parameters::getSchemeFileName()
{
    return scheme_file;
}


void Parameters::readObstacles(ifstream& in)
{
    string tmp="";
    unsigned num_obstacles = 0;

    while( !(str_dist(tmp,"</obstacle>") <= 2)){
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"<cylinders_list>") <= 2){
            readCylinderList(in);
            num_obstacles++;
        }
        if(str_dist(tmp,"<axons_list>") <= 2){
            readAxonList(in);
            num_obstacles++;
        }
        if((str_dist(tmp,"<neurons_list>") <= 2)){
            readNeuronList(in);
            num_obstacles++;
        }
        if(str_dist(tmp,"oriented_cylinders_list") <= 2){
            string path;
            in >> path;
            cylinders_files.push_back(path);
            num_obstacles++;
        }
        if(str_dist(tmp,"<ply>") <= 1){
            readPLY(in);
            num_obstacles++;
        }

        if(str_dist(tmp,"<new_packing>") <=1){
            readPackingParams(in);
            num_obstacles++;
        }

        if(str_dist(tmp,"<spheres_list>") <= 2){
            readSphereList(in);
            num_obstacles++;   
        }        

    }
    if(num_obstacles ==0){
        SimErrno::warning("<obstacle> tag initialized, but no valid obstacle tag found",cout);
    }
}
void Parameters::readVoxels(ifstream& in)
{
    string tmp="";
    double x,y,z;

    in >> x >> y >> z;
    min_limits = {x,y,z};
    in >> x >> y >> z;
    max_limits = {x,y,z};

    pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(min_limits,max_limits);
    voxels_list.push_back(voxel_);

    in >> tmp;
}

int Parameters::str_dist(string s, string t)
{
    ulong len_s = s.length();
    ulong len_t = t.length();

    /* base case: empty strings */
    if (len_s == 0) return int(len_t);
    if (len_t == 0) return int(len_s);

    if(len_s == 1 && len_t ==1)
        return s[0] != t[0];

    Eigen::MatrixXd costos(len_s,len_t);

    for(unsigned i = 0 ; i < s.size(); i++){
        for (unsigned j = 0 ; j < t.size(); j++){
            costos(i,j) = 0;
            costos(0,j) = j;
        }
        costos(i,0) = i;
    }

    int cost;

    for(unsigned i = 1 ; i < s.size(); i++){
        for (unsigned j = 1 ; j < t.size(); j++){
            /* test if last characters of the strings match */
            if (s[i] == t[j])
                cost = 0;
            else
                cost = 1;

            /* return minimum of delete char from s, delete char from t, and delete char from both */
            costos(i,j) =  min(min( costos(i-1,j) + 1,
                                    costos(i,j-1) + 1),
                               costos(i-1,j-1) + cost);
        }
    }

    return costos(s.length()-1,t.length()-1);
}

void Parameters::readInfoGatheringParams(ifstream& in)
{
    string tmp="";
    while(str_dist(tmp,"</log>"))
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"<positions>") <= 3)
        {
            in >> tmp;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

            while(str_dist(tmp,"</positions>"))
            {
                if(tmp.compare("t") == 0){
                    record_pos_times.push_back(num_steps);
                }
                else{
                    unsigned long time = stoul(tmp);
                    record_pos_times.push_back(unsigned(time));
                }

                in >> tmp;
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            }
        }
        else if(str_dist(tmp,"<phase>") <= 2)
        {
            in >> tmp;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

            while(str_dist(tmp,"</phase>"))
            {
                if(tmp.compare("t") == 0){
                    record_phase_times.push_back(num_steps);
                }
                else{
                    unsigned time = unsigned(stoul(tmp));
                    record_phase_times.push_back(time);
                }
                in >> tmp;
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            }
        }
        else if(str_dist(tmp,"<propagator>")<= 2){
            in >> tmp;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

            while(str_dist(tmp,"</propagator>"))
            {
                if(str_dist(tmp,"directions")<2){
                    std::string dir_path;
                    in >> dir_path;
                    cout << dir_path << endl; 
                    readPropagatorDirections(dir_path);
                    // Read directions
                }
                else if(tmp.compare("t") == 0){
                    unsigned step_idx;
                    in >> step_idx;
                    //this->record_prop_times.push_back(num_steps);
                    this->record_prop_times.push_back(step_idx);                    
                }
                else{
                    unsigned time = unsigned(stoul(tmp));
                    this->record_prop_times.push_back(time);
                }
                in >> tmp;
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

                std::sort(this->record_prop_times.begin(), this->record_prop_times.end());

                this->record_prop_times.erase(unique( this->record_prop_times.begin(), this->record_prop_times.end() ), this->record_prop_times.end() );
            }
            this->log_propagator = true;
        }
    }
}

void Parameters::readHexagonalParams(ifstream &in)
{
    hex_packing = true;

    string tmp="";

    while(str_dist(tmp,"</hex_packing>"))
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"radius") <= 1){
            in >> hex_packing_radius;
        }
        if(str_dist(tmp,"separation") <= 1){
            in >> hex_packing_separation;
        }
    }
}


void Parameters::readPackingParams(ifstream &in)
{
    string tmp="";

    while(str_dist(tmp,"</new_packing>") > 0)
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"<gamma_packing>") <= 1){
            readGammaParams(in);

        }
        if(str_dist(tmp,"<uniform_packing>") <= 1){
            readUniformParams(in);

        }
        if(str_dist(tmp,"<gaussian_packing>") <= 1){
            readGaussianParams(in);

        }
        if(str_dist(tmp,"<neuron_packing>") <= 1){
            readNeuronParams(in);
        }
        if(str_dist(tmp,"<hex_packing>") <=1){
            readHexagonalParams(in);
        }

        if(str_dist(tmp,"obstacle_type") <= 1){
            readObstacleType(in);
        }
        else if(str_dist(tmp,"</new_packing>") ==0){
            break;
        }

        tmp = "";
    }
}

void Parameters::readObstacleType(ifstream &in)
{
    string tmp;
    in >> tmp;
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
    if(str_dist(tmp,"sphere") <= 1){
        packing_s = true;
    }
    else if(str_dist(tmp,"cylinder") <= 1){
        packing_cyl = true;
    }
}

void Parameters::readGammaParams(ifstream &in)
{

    gamma_packing = true;

    string tmp="";

    while(str_dist(tmp,"</gamma_packing>") > 0)
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
        if(str_dist(tmp,"output_conf") <= 1){
            string tst;
            in >> tst;
            in >> packing_output_conf;
        }
        else if(str_dist(tmp,"alpha") <= 1 or str_dist(tmp,"shape") <= 1){
            string tmp2 ;
            getline(in, tmp2);
            istringstream alphas(tmp2);
            double alpha_;
            while(alphas >> alpha_){gamma_packing_alpha.push_back(alpha_);}
        }
        else if(str_dist(tmp,"beta") <= 1 or str_dist(tmp,"scale") <= 1){
            string tmp2 ;
            getline(in, tmp2);
            istringstream betas(tmp2);
            double beta_;
            while(betas >> beta_){gamma_packing_beta.push_back(beta_);}
        }
        else if(str_dist(tmp,"num_obstacles") <= 1){
            string tmp2 ;
            getline(in, tmp2);
            istringstream nbsphs(tmp2);
            unsigned nbsph_;
            while(nbsphs >> nbsph_){packing_num_obstacles.push_back(nbsph_);}
        }
        else if(str_dist(tmp,"icvf") <= 1){
            in >> packing_icvf;
        }
        else if(str_dist(tmp,"") == 0){
            in.clear();
            //in.ignore();
        }
        else if(str_dist(tmp,"</gamma_packing>") ==0){
            break;
        }

        tmp = "";
    }
}


void Parameters::readUniformParams(ifstream &in)
{

    uniform_packing = true;

    string tmp="";

    while(str_dist(tmp,"</uniform_packing>") > 0)
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
        if(str_dist(tmp,"output_conf") <= 1){
            string tst;
            in >> tst;
            in >> packing_output_conf;
        }
        else if(str_dist(tmp,"num_obstacles") <= 1){
    
            string line;
            getline(in, line);
            istringstream iss(line);
            unsigned subs;
            while (iss >> subs){packing_num_obstacles.push_back(subs);};

        }
        else if(str_dist(tmp,"radii") <= 1){
            
            string line="";
            getline(in, line);
            istringstream iss(line);
            double subs;
            while (iss >> subs){uniform_packing_radii.push_back(subs);};

        }        
        else if(str_dist(tmp,"icvf") <= 1){
            in >> packing_icvf;
        }
        else if(str_dist(tmp,"") == 0){
            in.clear();
            //in.ignore();
        }
        else if(str_dist(tmp,"</uniform_packing>") ==0){
            break;
        }

        tmp = "";
    }
}

void Parameters::readGaussianParams(ifstream &in)
{

    gaussian_packing = true;

    string tmp="";

    while(str_dist(tmp,"</gaussian_packing>") > 0)
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
        if(str_dist(tmp,"output_conf") <= 1){
            string tst;
            in >> tst;
            in >> packing_output_conf;
        }
        else if(str_dist(tmp,"num_obstacles") <= 1){
    
            string line;
            getline(in, line);
            istringstream iss(line);
            unsigned subs;
            while (iss >> subs){packing_num_obstacles.push_back(subs);};

        }
        else if(str_dist(tmp,"mean") <= 1){
            
            string line="";
            getline(in, line);
            istringstream iss(line);
            double subs;
            while (iss >> subs){gaussian_packing_mean.push_back(subs);};

        } 
        else if(str_dist(tmp,"std") <= 1){
            
            string line="";
            getline(in, line);
            istringstream iss(line);
            double subs;
            while (iss >> subs){gaussian_packing_std.push_back(subs);};

        }        
        else if(str_dist(tmp,"icvf") <= 1){
            in >> packing_icvf;
        }
        else if(str_dist(tmp,"") == 0){
            in.clear();
            //in.ignore();
        }
        else if(str_dist(tmp,"</gaussian_packing>") ==0){
            break;
        }

        tmp = "";
    }
}

void Parameters::readNeuronParams(ifstream &in)
{

    neuron_packing = true;

    string tmp="";

    while(str_dist(tmp,"</neuron_packing>") > 0)
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"num_obstacles") <= 1){
    
            string line;
            getline(in, line);
            istringstream iss(line);
            unsigned subs;
            while (iss >> subs){num_neurons = subs;};

        }      
        else if(str_dist(tmp,"icvf") <= 1){
            in >> packing_icvf;
        }
        else if(str_dist(tmp,"") == 0){
            in.clear();
            //in.ignore();
        }
        else if(str_dist(tmp,"permeability") <= 2){
        
            in >> tmp;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            if(str_dist(tmp,"global") <= 1){
                // One permeability for all obstacles
                in >> obstacle_permeability;
            }
            if(str_dist(tmp,"local") <= 1){
                // One permeability per obstacles
                string path;
                in >> path;
                neuron_permeability_files.push_back(path);
            }
        }
        else if(str_dist(tmp,"</neuron_packing>") == 0){
            break;
        }

        tmp = "";
    }
}

void Parameters::readSubdivisionFile()
{
    ifstream in(subdivisions_file);

    Subdivision tmp;

    while( in >> tmp.min_limits[0]){
        in >> tmp.min_limits[1];
        in >> tmp.min_limits[2];
        in >> tmp.max_limits[0];
        in >> tmp.max_limits[1];
        in >> tmp.max_limits[2];
        subdivisions.push_back(tmp);
    }
}


void Parameters::addSubdivisions()
{

    if( (number_subdivisions > 0) && (voxels_list.size() <=0) ){
        SimErrno::error("subdivisions_number parameter passed without a defined voxel.",cout);
        assert(0);
    }
    float gap[3];
    for(uint i = 0; i < 3; i++){
        gap[i] =  float(this->voxels_list[0].second[i] - this->voxels_list[0].first[i])/float(number_subdivisions);
    }

    for(uint x_ind = 0 ; x_ind < number_subdivisions; x_ind++){
        for(uint y_ind = 0 ; y_ind < number_subdivisions; y_ind++){
            for(uint z_ind = 0; z_ind < number_subdivisions; z_ind++){

                Subdivision tmp;

                //x
                tmp.min_limits[0] = float(voxels_list[0].first[0]  + x_ind*gap[0]);
                //y
                tmp.min_limits[1] = float(voxels_list[0].first[1]  + y_ind*gap[1]);
                //z
                tmp.min_limits[2] = float(voxels_list[0].first[2]  + z_ind*gap[2]);

                //x
                tmp.max_limits[0] = tmp.min_limits[0] + gap[0];
                //y
                tmp.max_limits[1] = tmp.min_limits[1] + gap[1];
                //z
                tmp.max_limits[2] = tmp.min_limits[2] + gap[2];

                subdivisions.push_back(tmp);
            }
        }
    }
}

void Parameters::readPropagatorDirections(string dir_path)
{
    ifstream in(dir_path);

    if(in.fail()){

        cout << " WHY!!! " << endl;
        //File does not exist code here
    }

    Eigen::Vector3f direction;


    while( in >> direction[0]){
        in >> direction[1];
        in >> direction[2];

        if(direction.norm() >0)
            this->prop_dirs.push_back(direction.normalized());
    }

    cout << direction << endl;

    in.close();
}


void Parameters::readPLY(ifstream& in)
{
    string path;
    in >> path;
    PLY_files.push_back(path);

    // Default scale is 1 and permeability 0 - Need one parameter per PLY obstacle for "mcsimulation" class

    string tmp      = "";
    double perm_    = 0.0;
    double scale_   = 1.0;

    while(!(str_dist(tmp,"</ply>") <= 2)){
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
        // One permeability per PLY file
        //in >> obstacle_permeability;
        if((str_dist(tmp,"permeability") <= 1)){in >> perm_;}
        if(str_dist(tmp,"ply_scale") <= 1){in >> scale_;}
    }
    PLY_scales.push_back(scale_);
    PLY_permeability.push_back(perm_);
         
}

void Parameters::readSphereList(ifstream& in)
{
    string path;
    in >> path;
    spheres_files.push_back(path);

    string tmp="";
    while(!(str_dist(tmp,"</spheres_list>") <= 2)){
        in >> tmp;
        if(!(str_dist(tmp,"permeability") <= 2)){
        
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            if(str_dist(tmp,"global") <= 1){
                // One permeability for all obstacles
                in >> obstacle_permeability;
            }
            if(str_dist(tmp,"local") <= 1){
                // One permeability per obstacles
                string path;
                in >> path;
                sphere_permeability_files.push_back(path);
            }
        }
    } 
}

void Parameters::readCylinderList(ifstream& in)
{
    string path;
    in >> path;
    cylinders_files.push_back(path);

    string tmp="";
    while(!(str_dist(tmp,"</cylinder_list>") <= 2)){
        in >> tmp;    
        if(!(str_dist(tmp,"permeability") <= 2)){
                    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            if(str_dist(tmp,"global") <= 1){
                // One permeability for all obstacles
                in >> obstacle_permeability;
            }
            if(str_dist(tmp,"local") <= 1){
                // One permeability per obstacles
                string path;
                in >> path;
                cylinder_permeability_files.push_back(path);
            }
        }
    } 
}

void Parameters::readAxonList(ifstream& in)
{
    string path;
    in >> path;
    axons_files.push_back(path);

    string tmp="";
    while(!(str_dist(tmp,"</axons_list>") <= 2)){
        in >> tmp;    
        if(!(str_dist(tmp,"permeability") <= 2)){
                    std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            if(str_dist(tmp,"global") <= 1){
                // One permeability for all obstacles
                in >> obstacle_permeability;
            }
            if(str_dist(tmp,"local") <= 1){
                // One permeability per obstacles
                string path;
                in >> path;
                axon_permeability_files.push_back(path);
            }
        }
    } 
}

void Parameters::readNeuronList(ifstream& in)
{
    string path;
    in >> path;
    neurons_files.push_back(path);

    string tmp="";

    while(!(str_dist(tmp,"</neurons_list>") <= 2))
    {
        in >> tmp;  

        if(!(str_dist(tmp,"permeability") <= 2))
        {
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            if(str_dist(tmp,"global") <= 1)
            {
                // One permeability for all obstacles
                in >> obstacle_permeability;
            }
            if(str_dist(tmp,"local") <= 1)
            {
                // One permeability per obstacles
                string path;
                in >> path;
                neuron_permeability_files.push_back(path);
            }
        }
    } 
}