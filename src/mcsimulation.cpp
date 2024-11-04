#include "mcsimulation.h"
#include <Eigen/Dense>
#include "simerrno.h"
#include "pgsesequence.h"
#include "gradientwaveform.h"
#include <iostream>

int MCSimulation::count =0;

MCSimulation::MCSimulation()
{
    dynamicsEngine = NULL;
    dataSynth = NULL;
    id = count;
    count++;
}

/*DEPRECATED*/
MCSimulation::MCSimulation(std::string config_file)
{
    dynamicsEngine = NULL;
    dataSynth      = NULL;

    params.readSchemeFile(config_file);
    dynamicsEngine = new DynamicsSimulation(params);

    if(params.scheme_file.length() > 2){
        scheme.readSchemeFile(params.scheme_file,params.scale_from_stu);
    }

    if(scheme.type == "PGSE"){
        dataSynth = new PGSESequence(scheme);
        dataSynth->setNumberOfSteps(dynamicsEngine->params.num_steps);

        if(params.subdivision_flag){
            dataSynth->subdivision_flag = true;
            dataSynth->subdivisions = params.subdivisions;
            dataSynth->initializeSubdivisionSignals();
        }
    }

    dynamicsEngine->id = count;
    id = count;
    count++;
}

MCSimulation::MCSimulation(Parameters& params_)
{
    dynamicsEngine = NULL;
    dataSynth = NULL;

    params = params_;

    dynamicsEngine = new DynamicsSimulation(params);
  

    if(params.scheme_file.length() > 2){
        scheme.readSchemeFile(params.scheme_file,params.scale_from_stu);
    }

    if(scheme.type == "PGSE"){
        dataSynth = new PGSESequence(scheme);
    }
    if(scheme.type == "WAVEFORM"){
        dataSynth = new GradientWaveform(scheme);
    }

    if (dataSynth){
        dataSynth->setNumberOfSteps(dynamicsEngine->params.num_steps);
    }
    
    if(params.subdivision_flag){
        dataSynth->subdivision_flag = true;
        dataSynth->subdivisions = params.subdivisions;
        dataSynth->initializeSubdivisionSignals();
    }

    dynamicsEngine->id = count;
    id = count;
    count++;
}


void MCSimulation::startSimulation()
{

    iniObstacles();
    // update number of walkers
    dynamicsEngine->params = params;
    //cout << " number of walkers : " << params.num_walkers << endl;

    if(dataSynth != NULL){
        dynamicsEngine->startSimulation(dataSynth);
    }
    else{
        dynamicsEngine->startSimulation();
    }

}

double MCSimulation::getExpectedFreeeDecay(unsigned i)
{
    if(dataSynth){
        double b = dataSynth->getbValue(i);
        return exp(-b*params.diffusivity_extra);
    }

    return -1;
}


void MCSimulation::iniObstacles()
{
    addCylindersObstaclesFromFiles();

    addAxonsObstaclesFromFiles();

    addNeuronsObstaclesFromFiles();

    addPLYObstaclesFromFiles();

    addVoxels();

    addCylindersConfigurations();
    //Used only if there's a voxel (deprecated)
    //addExtraObstacles();

    addSpheresObstaclesFromFiles();

}


//* Auxiliare method to split words in a line using the spaces*//
template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}



bool withinBounds(Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, Eigen::Vector3d pos, double distance)
{
    bool within;
    for (int i = 0; i < 3; i++) // check for all dimensions
    {
        if ((pos[i] < max_limits[i] + distance) && (pos[i] > min_limits[i] - distance))
        {
            within = true;
        }
        else
        {
            within = false;
            break;
        }
    }
    return within;
}


bool withinAreaBounds(Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, Eigen::Vector3d pos, double distance)
{
    bool within;
    for (int i = 0; i < 2; i++) // check for all dimensions
    {
        if ((pos[i] < max_limits[i] + distance) && (pos[i] > min_limits[i] - distance))
        {
            within = true;
        }
        else
        {
            within = false;
            break;
        }
    }
    return within;
}


double computeAreaICVF(Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector <Cylinder> cylinders_)
{
    if (cylinders_.size() == 0)
        return 0;
    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]); // total area
    double AreaC = 0;

    for (uint i = 0; i < cylinders_.size(); i++) // for all axons
    {

        if (withinAreaBounds(min_limits, max_limits, cylinders_[i].P, cylinders_[i].radius))
        {
   
            AreaC +=  M_PI * cylinders_[i].radius * cylinders_[i].radius;
        }
        else if (withinAreaBounds(min_limits, max_limits, cylinders_[i].P, 0))
        {
   
            AreaC +=  M_PI * cylinders_[i].radius * cylinders_[i].radius/2;
        }
    }
    return AreaC / AreaV; // ( total axons volume / total volume )
}

double computeICVF(Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector <Axon> axons)
{
    if (axons.size() == 0)
        return 0;
    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]) * (max_limits[2] - min_limits[2]); // total volume
    double AreaC = 0;

    for (uint i = 0; i < axons.size(); i++) // for all axons
    {

        if (axons[i].spheres.size() > 1)
        {
            for (uint j = 1; j < axons[i].spheres.size(); j++)
            {
                double l = (axons[i].spheres[j - 1].P - axons[i].spheres[j].P).norm(); // distance between centers
                double mean_r = (axons[i].spheres[j - 1].radius + axons[i].spheres[j].radius) / 2;

                if (withinBounds(min_limits, max_limits,axons[i].spheres[j].P, axons[i].spheres[j].radius) && withinBounds(min_limits, max_limits,axons[i].spheres[j-1].P, axons[i].spheres[j-1].radius))
                {
                    AreaC += l * M_PI * mean_r * mean_r;
                }
                else if (withinBounds(min_limits, max_limits,axons[i].spheres[j].P, 0) && withinBounds(min_limits, max_limits,axons[i].spheres[j-1].P, 0))
                {
                    AreaC += l * M_PI * mean_r * mean_r/2;
                }
            }
        }
    }
    return AreaC / AreaV; // ( total axons volume / total volume )
}


int MCSimulation::str_dist(string s, string t)
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


void MCSimulation::addAxonsObstaclesFromFiles()
{
    
    for(unsigned i = 0; i < params.axons_files.size(); i++){
        cout << "Adding Axons" << endl;


        std::ifstream in(params.axons_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}

            std::vector<std::string> jkr = split(line,' ');
            if (jkr.size() != 10){
                //std::cout << "\033[1;33m[Warning]\033[0m Cylinder orientation was set towards the Z direction by default" << std::endl;
            }
            break;
        }
        in.close();

        // Permeability file - if any
        double perm_; 

        std::ifstream in_perm;
        if(params.axon_permeability_files.size() >0){
            in_perm.open(params.axon_permeability_files[i]);
        }

        // Diffusion coefficients
        double diff_i; 
        double diff_e;

        in.open(params.axons_files[i]);
        double x,y,z,rout, rin, p, r;
        double ax_id_, sph_id_, branch_id_;
        int ax_id, sph_id, branch_id;
        int last_ax_id = -1;
        std::string type_object, last_type ="";
        std::string header;

        std::vector<Sphere> spheres_out ;
        std::vector<Sphere> spheres_in ;
        Sphere sphere_out;
        Sphere sphere_in;

        int line_num = 0;

        int header_size = 10;

        for(unsigned j = 0; j < header_size; j++){  
            in >>header;
            //cout << "header :" << header << endl;
        } 

        while (in >>ax_id_ >> sph_id_ >> branch_id_ >> type_object >> x >> y >> z >> rin >> rout >> p){
         
            // convert um to m
            x = x/1000.0;
            y = y/1000.0;
            z = z/1000.0;
            rout = rout/1000.0;
            rin = rin/1000.0;

            sph_id = int(sph_id_);
            ax_id = int(ax_id_);
            branch_id = int(branch_id_);
            //cout << "x :" << x << " y :" << y << " z :" << z << " rin :" << rin << " rout :" << rout << endl;

            // if the new line is from a different axon
            if (line_num !=0 and last_ax_id != ax_id and str_dist(last_type,"axon") <= 1){
                // create the axon with id : last_ax_id
                Axon ax (last_ax_id, {0.0,0.0,0.0}, {0.0,0.0,0.0}, rout);
                Axon ax_in (last_ax_id, {0.0,0.0,0.0}, {0.0,0.0,0.0}, rin);

                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  
                
                //cout << "perm_ :"   << perm_ << endl;
                for (unsigned i = 0; i < spheres_out.size(); i++){
                    spheres_out[i].setPercolation(perm_);
                    spheres_in[i].setPercolation(perm_);
                    // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                    diff_i = params.diffusivity_intra; 
                    diff_e = params.diffusivity_extra;
                    spheres_out[i].setDiffusion(diff_i, diff_e);
                    spheres_in[i].setDiffusion(diff_i, diff_e);
                }
                ax.setDiffusion(diff_i, diff_e);
                ax.setPercolation(perm_);
                ax.set_spheres(spheres_out);
                dynamicsEngine->axons_list.push_back(ax);

                ax_in.setDiffusion(diff_i, diff_e);
                ax_in.setPercolation(perm_);
                ax_in.set_spheres(spheres_in);
                dynamicsEngine->inner_axons_list.push_back(ax_in);
                spheres_out.clear();
                spheres_in.clear();

            }
            sphere_out = Sphere(sph_id, ax_id, Eigen::Vector3d(x,y,z), rout, 0);
            sphere_in = Sphere(sph_id, ax_id, Eigen::Vector3d(x,y,z), rin, 0);
            spheres_out.push_back(sphere_out);
            spheres_in.push_back(sphere_in);
            last_ax_id = ax_id;
            last_type = type_object;
            line_num += 1;
            
                
        }
        
        if (str_dist(last_type,"axon") <= 1) {
            // add last sphere
            Axon ax (last_ax_id, {0.0,0.0,0.0}, {0.0,0.0,0.0}, rout);
            Axon ax_in (last_ax_id, {0.0,0.0,0.0}, {0.0,0.0,0.0}, rin);

            for (unsigned i = 0; i < spheres_out.size(); i++){
                spheres_out[i].setPercolation(perm_);
                spheres_in[i].setPercolation(perm_);
                // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                diff_i = params.diffusivity_intra; 
                diff_e = params.diffusivity_extra;
                spheres_out[i].setDiffusion(diff_i, diff_e);
                spheres_in[i].setDiffusion(diff_i, diff_e);
            }
            ax.setDiffusion(diff_i, diff_e);
            ax.setPercolation(perm_);
            ax.set_spheres(spheres_out);
            dynamicsEngine->axons_list.push_back(ax);

            ax_in.setDiffusion(diff_i, diff_e);
            ax_in.setPercolation(perm_);
            ax_in.set_spheres(spheres_in);
            dynamicsEngine->inner_axons_list.push_back(ax_in);
            spheres_out.clear();
            spheres_in.clear();


        }

        //double max_limits, min_limits;
        // z of last sphere of first axon
        //Axon first_axon = dynamicsEngine->axons_list[0];
        //cout << "first_axon.spheres.size() : " << first_axon.spheres.size() << endl;
        //max_limits = first_axon.spheres[first_axon.spheres.size()-1].P[2]; 
        // z of first sphere of first axon
        //min_limits = first_axon.spheres[0].P[2]; 
        //max_limits =0.05; 
        //cout << "min_limits : " << min_limits << endl;
        //cout << "max_limits : " << max_limits << endl;

        // voxel
        //params.max_limits = Eigen::Vector3d(max_limits,max_limits,max_limits);
        //params.min_limits = Eigen::Vector3d(min_limits,min_limits,min_limits);
        //double volume = (params.max_limits[0]-params.min_limits[0])*(params.max_limits[1]-params.min_limits[1])*(params.max_limits[2]-params.min_limits[2]);
        //pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(min_limits,max_limits);
        //params.voxels_list.push_back(voxel_);
    
        //cout << " Volume :" << volume << endl;
        // set icvf
        //double icvf = computeICVF(params.min_limits, params.max_limits, dynamicsEngine->axons_list);
        // set number of particles
        //if (params.concentration != 0){
        //    if (params.ini_walker_flag == "intra"){
        //        int num_walkers = params.concentration*volume*icvf/params.num_proc;
        //        params.setNumWalkers(num_walkers);
        //    }
        //    else if (params.ini_walker_flag == "extra"){
        //        int num_walkers = params.concentration*volume*(1.0-icvf)/params.num_proc;
        //        params.setNumWalkers(num_walkers);
        //    }
        //    else{
        //        int num_walkers = params.concentration*volume/params.num_proc;
        //        params.setNumWalkers(num_walkers);
        //    } 
        //}
        cout << "params.ini_walker_flag :" << params.ini_walker_flag << endl;
        
        //cout << " ICVF :" << icvf<< endl;
        cout << " Number of particles :" << params.num_walkers << endl;
        cout << "Number of axons :" << dynamicsEngine->axons_list.size() << endl;

        in.close();
        
    }
}




// void MCSimulation::addGlialsObstaclesFromFiles()
// {
    

//     for(unsigned i = 0; i < params.glials_files.size(); i++){
//         cout << "Adding Glials" << endl;


//         std::ifstream in(params.glials_files[i]);

//         if(!in){
//             return;
//         }

//         bool first=true;
//         for( std::string line; getline( in, line ); )
//         {
//             if(first) {first-=1;continue;}

//             std::vector<std::string> jkr = split(line,' ');
//             if (jkr.size() != 10){
//                 //std::cout << "\033[1;33m[Warning]\033[0m Cylinder orientation was set towards the Z direction by default" << std::endl;
//             }
//             break;
//         }
//         in.close();

//         // Permeability file - if any
//         double perm_; 

//         std::ifstream in_perm;
//         if(params.glial_permeability_files.size() >0){
//             in_perm.open(params.glial_permeability_files[i]);
//         }

//         // Diffusion coefficients
//         double diff_i; 
//         double diff_e;

//         in.open(params.glials_files[i]);
//         double x,y,z,rout, rin, p, r;
//         int ax_id, sph_id, branch_id;
//         std::string type_object;
//         std::string header;

//         std::vector<Sphere> spheres_ ;
//         Sphere sphere_;
//         int line_num = 0;

//         int header_size = 10;

//         for(unsigned j = 0; j < header_size; j++){  
//             in >>header;
//             cout << "header :" << header << endl;
//         } 
//         std::vector <Sphere> processes_ = std::vector<Sphere>();
//         Glial glial_cell;
//         diff_i = params.diffusivity_intra; 
//         diff_e = params.diffusivity_extra;
//         perm_ = params.obstacle_permeability;


//         while (in >>ax_id >> sph_id >> branch_id>> type_object >> x >> y >> z >> rin >> rout >> p){

//             x = x/1000.0;
//             y = y/1000.0;
//             z = z/1000.0;
//             r = rout/1000.0;
            
//             if (str_dist(type_object,"glial") == 0){
//                 //cout << endl;
//                 //cout << "processes_.size() :" << processes_.size() << endl;
//                 if (processes_.size() != 0){
//                     // create the glial with id : last_ax_id
//                     glial_cell.setDiffusion(diff_i, diff_e);
//                     glial_cell.setPercolation(perm_);
//                     //cout << "processes_.size() added to glial :" << processes_.size() << " line num :" << line_num << endl;
//                     glial_cell.set_spheres(processes_);
//                     processes_.clear();
//                     processes_ = std::vector<Sphere>();
//                     dynamicsEngine->glials_list.push_back(glial_cell);
//                 }
     
//                 Sphere soma = Sphere (sph_id, ax_id, Eigen::Vector3d(x,y,z), r, 1);
//                 soma.setDiffusion(diff_i, diff_e);
//                 soma.setPercolation(perm_);
//                 glial_cell = Glial(ax_id, soma);
//             }
//             if (str_dist(type_object,"glialRamification") == 0){
        
//                 Sphere process = Sphere (sph_id, ax_id, Eigen::Vector3d(x,y,z), r, 1);
//                 process.setDiffusion(diff_i, diff_e);
//                 process.setPercolation(perm_);
//                 processes_.push_back(process);
//                 //cout << "processes_.size() :" << processes_.size() << endl;
//             }
//             line_num += 1;

//         }
//         //last cell
//         glial_cell.setDiffusion(diff_i, diff_e);
//         glial_cell.setPercolation(perm_);
//         //cout << "processes_.size() added to glial :" << processes_.size() << " line num :" << line_num << endl;
//         glial_cell.set_spheres(processes_);
//         processes_.clear();
//         dynamicsEngine->glials_list.push_back(glial_cell);

//         cout << "Number of glials :" << dynamicsEngine->glials_list.size() << endl;
//         cout << dynamicsEngine->glials_list[0].processes.size() << endl;
//         cout << dynamicsEngine->glials_list[0].processes[0].P  << endl;

//         in.close();
//     }
// }

void MCSimulation::addNeuronsObstaclesFromFiles()
{
    

    for(unsigned i = 0; i < params.glials_files.size(); i++){
        std::ifstream in(params.glials_files[i]);

        if(!in){
            std::cout <<  "[ERROR] Unable to open:" << params.glials_files[i] << std::endl;
            return;
        }
        unsigned enum_ = 1;

        bool first = true;

        for( std::string line; getline( in, line ); )
        {
            if(first) {
                first  = false;
                enum_ += 1;
                continue;
                }
            if (enum_ == 2 || enum_ == 3 || enum_ == 4 || enum_ == 5 || enum_ == 6){
                enum_ += 1;
                continue;
            }

            std::vector<std::string> jkr = split(line,' ');
            if (jkr.size() != 4 && jkr.size() != 2 && jkr.size() != 7){
                std::cout << jkr.size() <<  " elements per line" << std::endl;
                std::cout << "wrong number of elements per line in file" << std::endl;
            }
            break;
        }
        in.close();
        // Permeability file - if any
        double perm_; 

        std::ifstream in_perm;
        if(params.glial_permeability_files.size() > 0)
            in_perm.open(params.glial_permeability_files[i]);
        

        // Local permeability - Different for each obstacle
        if(params.glial_permeability_files.size() > 0)
            in_perm >> perm_;
        // Global permeability - Same for all obstacle
        else
            perm_ = params.obstacle_permeability;
        
                
        // Diffusion coefficients
        double diff_i; 
        double diff_e;

        // cout << params.neurons_files[neurons_files_id] << endl;
        double x,y,z,r;
        double running_time, icvf;

        in.open(params.glials_files[i]);

        string part;
        std::string line;
        int sph_id = 0;

        Glial glial_cell;
        std::vector<std::vector <Sphere>> processes_;
        std::vector <Sphere> process;
        int neuron_id, dendrite_id = -1;
        int last_dendrite_id = 0;

        while (getline( in, line )) {
            string type;
    
            if(line.size() > 0)
            {
                vector<string> jkr = split(line,' ');
                // for(int i=0; i < jkr.size(); i++)
                //     cout << jkr[i] << " ";
                // cout << endl;
                // Local permeability - Different for each obstacle
                if(in_perm)
                    in_perm >> perm_;
                // Global permeability - Same for all obstacle
                else
                    perm_ = params.obstacle_permeability; 

                // Process neuron
                if (jkr[0] == "Neuron") 
                {
                    getline( in, line );
                    vector<string> jkr = split(line,' ');
                    neuron_id   = stod(jkr[0]);
                    dendrite_id = stod(jkr[1]);
                    x = stod(jkr[3]) / 1000.0;
                    y = stod(jkr[4]) / 1000.0;
                    z = stod(jkr[5]) / 1000.0;
                    r = stod(jkr[6]) / 1000.0;

                    if(r > 0)
                    {
                        Sphere soma = Sphere (0, processes_.size(), Eigen::Vector3d(x,y,z), r, 1);
                        soma.setDiffusion(diff_i, diff_e);
                        soma.setPercolation(perm_);
                        glial_cell = Glial(dynamicsEngine->glials_list.size(), soma);
                    }
                    else
                    {
                        Sphere soma = Sphere (0, processes_.size(), Eigen::Vector3d(x,y,z), 0, 1);
                        soma.setDiffusion(diff_i, diff_e);
                        soma.setPercolation(perm_);
                        glial_cell = Glial(dynamicsEngine->glials_list.size(), soma);
                    }
                    
                } 
                // end neuron
                else if (jkr[0] == "end") 
                {
                    glial_cell.set_spheres(process);
                    // create the glial with id : last_ax_id
                    glial_cell.setDiffusion(diff_i, diff_e);
                    glial_cell.setPercolation(perm_);
                    glial_cell.set_box_neuron();
                    dynamicsEngine->glials_list.push_back(glial_cell);
                    dendrite_id = last_dendrite_id = -1;
                    processes_.clear();
                    process.clear();
                    glial_cell.Box.clear();
                    glial_cell.Box_branch.clear();
                } // end neuron
                // Process dendrite
                else if (jkr.size() > 3)  
                {
                    neuron_id   = stod(jkr[0]);
                    dendrite_id = stod(jkr[1]);
                    x = stod(jkr[3]) / 1000.0;
                    y = stod(jkr[4]) / 1000.0;
                    z = stod(jkr[5]) / 1000.0;
                    r = stod(jkr[6]) / 1000.0;
                    // cout << x << " " << y << " " << z << " " << r << endl;
                    
                    // New dendrite
                    if(dendrite_id != last_dendrite_id)
                    {
                        sph_id = 0;
                        if(last_dendrite_id != -1)
                        {
                            glial_cell.set_spheres(process);
                            process.clear();
                        }
                        last_dendrite_id = dendrite_id;
                    }
                    Sphere sph = Sphere(sph_id, dendrite_id, Eigen::Vector3d(x,y,z), r, 1);
                    sph.setDiffusion(diff_i, diff_e);
                    sph.setPercolation(perm_);
                    process.push_back(sph);
                    sph_id++;
                }    
            } // if( line.size() > 0)
        
        }//while( getline(in, line) )
            
        params.gamma_icvf = icvf;

        in.close();

        double volume     = (params.max_limits[0] - params.min_limits[0]) * (params.max_limits[1] - params.min_limits[1]) * (params.max_limits[2] - params.min_limits[2]);
        // double icvf_calculated = computeICVF(params.min_limits, params.max_limits, dynamicsEngine->neurons_list);
        cout << "Number of glials :" << dynamicsEngine->glials_list.size() << endl;
        // cout << dynamicsEngine->glials_list[0].processes.size() << endl;
        // cout << dynamicsEngine->glials_list[0].processes[0].P  << endl;
    }
}

void MCSimulation::addCylindersObstaclesFromFiles()
{

   
    for(unsigned i = 0; i < params.cylinders_files.size(); i++){
        cout << "Adding Cylinders" << endl;


        std::ifstream in(params.cylinders_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}

            std::vector<std::string> jkr = split(line,' ');
            if (jkr.size() != 10){
                std::cout << "\033[1;33m[Warning]\033[0m Cylinder file does not have 10 elements per line" << std::endl;
            }
            break;
        }
        in.close();

        // Permeability file - if any
        double perm_; 

        std::ifstream in_perm;
        if(params.cylinder_permeability_files.size() >0){
            in_perm.open(params.cylinder_permeability_files[i]);
        }

        // Diffusion coefficients
        double diff_i; 
        double diff_e;

        in.open(params.cylinders_files[i]);
        double x,y,z,rout, rin, p, r, last_z;
        double ax_id, sph_id, branch_id;
        int last_ax_id = -1;
        std::string type_object, last_type ="";
        std::string header;

        int line_num = 0;

        int header_size = 10;

        for(unsigned j = 0; j < header_size; j++){  
            in >>header;
            //cout << "header :" << header << endl;
        } 

        diff_i = params.diffusivity_intra; 
        diff_e = params.diffusivity_extra;
        perm_ = params.obstacle_permeability;


        while (in >>ax_id >> sph_id >> branch_id>> type_object >> x >> y >> z >> rin >> rout >> p){
            //cout << "Ax_id :" << ax_id << endl;
            // convert to mm
            x = x/1000.0;
            y = y/1000.0;
            z = z/1000.0;
            rout = rout/1000.0;
            rin = rin/1000.0;

            // if the new line is from a different axon
            if (line_num !=0 and last_ax_id != ax_id and str_dist(last_type,"axon") <= 1){
                //cout << "ax_id :" << ax_id << endl;
                // create the axon with id : last_ax_id
                Cylinder cyl (last_ax_id, {x,y,0.0}, {x,y,last_z}, rout);
                Cylinder cyl_in (last_ax_id, {x,y,0.0}, {x,y,last_z}, rin);

                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  
                
                cyl.setDiffusion(diff_i, diff_e);
                cyl.setPercolation(perm_);
                cyl_in.setDiffusion(diff_i, diff_e);
                cyl_in.setPercolation(perm_);
                dynamicsEngine->cylinders_list.push_back(cyl);

                dynamicsEngine->inner_cylinders_list.push_back(cyl_in);

            }
            last_ax_id = ax_id;
            last_type = type_object;
            line_num += 1;
            last_z = z;
                
        }
        
        if (str_dist(last_type,"axon") <= 1) {
            // add last sphere
            Cylinder cyl (last_ax_id, {x,y, 0.0}, {x,y,last_z}, rout);
            Cylinder cyl_in (last_ax_id, {x,y, 0.0}, {x,y,last_z}, rin);

            cyl.setDiffusion(diff_i, diff_e);
            cyl.setPercolation(perm_);
            dynamicsEngine->cylinders_list.push_back(cyl);

            cyl_in.setDiffusion(diff_i, diff_e);
            cyl_in.setPercolation(perm_);
            dynamicsEngine->inner_cylinders_list.push_back(cyl_in);

        }

        cout << "params.ini_walker_flag :" << params.ini_walker_flag << endl;
        
        //cout << " ICVF :" << icvf<< endl;
        cout << " Number of particles :" << params.num_walkers << endl;
        cout << "Number of cylinders :" << dynamicsEngine->cylinders_list.size() << endl;

        in.close();
    }
}


void MCSimulation::addPLYObstaclesFromFiles()
{
    for(unsigned i = 0; i < params.PLY_files.size(); i++){

        PLYObstacle ply_(params.PLY_files[i],params.PLY_scales[i]);

        // Permeability - Kept outside initialization to be consistent with cylinders and spheres. Easily moved to ply constructor. 
        double perm_; 
        perm_ = params.PLY_permeability[i];
        ply_.setPercolation(perm_);

        // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
        double diff_i; 
        double diff_e;
        
        diff_i = params.diffusivity_intra; 
        diff_e = params.diffusivity_extra;
        ply_.setDiffusion(diff_i, diff_e);

        // Add PLY to list
        dynamicsEngine->plyObstacles_list.push_back(ply_);
    }
}

void MCSimulation::addVoxels()
{
    for(unsigned i = 0 ; i < params.voxels_list.size(); i++){
        dynamicsEngine->voxels_list.push_back(Voxel(params.voxels_list[i].first,params.voxels_list[i].second));
    }
}

void MCSimulation::addCylindersConfigurations()
{

    if(params.hex_packing){
        double rad = params.hex_packing_radius,sep = params.hex_packing_separation;

        // h = sqrt(3)/2 * sep
        double h = 0.866025404*sep;

        dynamicsEngine->cylinders_list.push_back(Cylinder(0,Eigen::Vector3d(0,0,0),Eigen::Vector3d(0,0,1.0),rad));
        dynamicsEngine->cylinders_list.push_back(Cylinder(0,Eigen::Vector3d(sep,0,0),Eigen::Vector3d(sep,0,1.0),rad));

        dynamicsEngine->cylinders_list.push_back(Cylinder(0,Eigen::Vector3d(0,2.0*h,0),Eigen::Vector3d(0,2.0*h,1.0),rad));
        dynamicsEngine->cylinders_list.push_back(Cylinder(0,Eigen::Vector3d(sep,2.0*h,0),Eigen::Vector3d(sep,2.0*h,1.0),rad));

        dynamicsEngine->cylinders_list.push_back(Cylinder(0,Eigen::Vector3d(0.5*sep,h,0),Eigen::Vector3d(0.5*sep,h,1.0),rad));

        // To avoid problems with the boundaries
        dynamicsEngine->cylinders_list.push_back(Cylinder(0,Eigen::Vector3d(-0.5*sep,h,0),Eigen::Vector3d(-0.5*sep,h,1.0),rad));
        dynamicsEngine->cylinders_list.push_back(Cylinder(0,Eigen::Vector3d(1.5*sep,h,0),Eigen::Vector3d(1.5*sep,h,1.0),rad));

        if(dynamicsEngine->voxels_list.size()>0)
            dynamicsEngine->voxels_list.clear();

        dynamicsEngine->voxels_list.push_back(Voxel(Eigen::Vector3d(0,0,0),Eigen::Vector3d(sep,2.0*h,2.0*h)));

    }
}

void MCSimulation::addSpheresObstaclesFromFiles()
{
    for(unsigned i = 0; i < params.spheres_files.size(); i++){

        std::ifstream in(params.spheres_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}
            break;
        }
        in.close();

        // Permeability file - if any
        double perm_; 

        std::ifstream in_perm;
        if(params.sphere_permeability_files.size() >0){
            in_perm.open(params.sphere_permeability_files[i]);
        }

        // Diffusion coefficients
        double diff_i; 
        double diff_e;
            
        in.open(params.spheres_files[i]);
        double x,y,z,r;
        double scale;
        in >> scale;

        while (in >> x >> y >> z >> r)
        {
            Sphere sph(0,0,Eigen::Vector3d(x,y,z),r,scale);

            // Local permeability - Different for each obstacle
            if(in_perm.is_open()){
                in_perm >> perm_;
            }
            // Global permeability - Same for all obstacle
            else{
                perm_ = params.obstacle_permeability;
            }            

            sph.setPercolation(perm_);

            // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
            diff_i = params.diffusivity_intra; 
            diff_e = params.diffusivity_extra;
            sph.setDiffusion(diff_i, diff_e);
            
            // Add sphere to list
            dynamicsEngine->spheres_list.push_back(sph);     
        }
        in.close();
        in_perm.close();
    }
}

bool cylinderIsCloseBoundery(Cylinder& cyl, Eigen::Vector3d min_limits,Eigen::Vector3d max_limits){

    //gap to the boundary
    double gap = 1e-6;
    //3 dimensional vector
    for (int i = 0 ; i < 3; i++)
        if( (cyl.P[i] - cyl.radius - gap < min_limits[i]) || (cyl.P[i] + cyl.radius + gap  > max_limits[i]) )
            return true;

    return false;
}

void MCSimulation::addExtraObstacles()
{
    if(dynamicsEngine->voxels_list.size() == 0)
        return;

    std::vector<Eigen::Vector3d> multipliers;

    Eigen::Vector3d gap = params.max_limits - params.min_limits;

    for(int i = -1  ;i <= 1; i++)
        for(int j = -1  ;j <= 1; j++)
            for(int k = -1 ;k <= 1; k++){
                Eigen::Vector3d jkr(i*gap[0],j*gap[1],k*gap[2]);
                multipliers.push_back(jkr);
            }


    unsigned long cylinders_num = dynamicsEngine->cylinders_list.size();

    for (unsigned c = 0; c < cylinders_num ;c++)
        for (unsigned i = 0 ; i < multipliers.size(); i++)
            if(multipliers[i][0]!=0.0 || multipliers[i][1]!=0.0 || multipliers[i][2]!=0.0)
            {
                Eigen::Vector3d P_ = dynamicsEngine->cylinders_list[c].P;
                Eigen::Vector3d Q_ = dynamicsEngine->cylinders_list[c].Q;
                P_[0]+= multipliers[i][0];P_[1]+=multipliers[i][1];P_[2]+=multipliers[i][2];
                Q_[0]+= multipliers[i][0];Q_[1]+=multipliers[i][1];Q_[2]+=multipliers[i][2];
                Cylinder tmp_cyl(0,P_,Q_,dynamicsEngine->cylinders_list[c].radius);


                //if the obstacle is close enough
                //if (cylinderIsCloseBoundery(tmp_cyl,params.min_limits,params.max_limits))
                    dynamicsEngine->cylinders_list.push_back(tmp_cyl);
            }

}


MCSimulation::~MCSimulation()
{
    if(dynamicsEngine != NULL)
        delete dynamicsEngine;

    if(dataSynth != NULL)
        delete dataSynth;
}


