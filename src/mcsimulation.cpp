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


void MCSimulation::addAxonsObstaclesFromFiles()
{

    for(unsigned i = 0; i < params.axons_files.size(); i++){


        std::ifstream in(params.axons_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}

            std::vector<std::string> jkr = split(line,' ');
            if (jkr.size() != 8){
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

        double x,y,z,r;
        int ax_id, sph_id, p;
        int last_ax_id = -1;
        std::string type_object;
        std::string header;

        std::vector<Sphere> spheres_ ;
        Sphere sphere_;
        int line_num = 0;

        int header_size = 8;

        for(unsigned j = 0; j < header_size; j++){  
            in >>header;
            //cout << "header :" << header << endl;
        } 

        while (in >>ax_id >> sph_id >> type_object >> x >> y >> z >> r >> p){

            // convert to mm
            x = x/1000.0;
            y = y/1000.0;
            z = z/1000.0;
            r = r/1000.0;

                
            // if the new line is from a different axon
            if (line_num !=0 and last_ax_id != ax_id){
                // create the axon with id : last_ax_id
                Axon ax (last_ax_id, {0.0,0.0,0.0}, {0.0,0.0,0.0}, r);

                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  
                
                for (unsigned i = 0; i < spheres_.size(); i++){
                    spheres_[i].setPercolation(perm_);
                    // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                    diff_i = params.diffusivity_intra; 
                    diff_e = params.diffusivity_extra;
                    spheres_[i].setDiffusion(diff_i, diff_e);
                }
                ax.set_spheres(spheres_);
                dynamicsEngine->axons_list.push_back(ax);
                spheres_.clear();

            }
            sphere_ = Sphere(sph_id, ax_id,Eigen::Vector3d(x,y,z), r);
            spheres_.push_back(sphere_);
            last_ax_id = ax_id;
            line_num += 1;
                
        }
        
        

        // add last sphere
        Axon ax (last_ax_id, {0.0,0.0,0.0}, {0.0,0.0,0.0}, r);
        for (unsigned i = 0; i < spheres_.size(); i++){
            spheres_[i].setPercolation(perm_);
            // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
            diff_i = params.diffusivity_intra; 
            diff_e = params.diffusivity_extra;
            spheres_[i].setDiffusion(diff_i, diff_e);
        }
        ax.set_spheres(spheres_);
        dynamicsEngine->axons_list.push_back(ax);

        double max_limits, min_limits;
        // z of last sphere of first axon
        Axon first_axon = dynamicsEngine->axons_list[0];
        max_limits = first_axon.spheres[first_axon.spheres.size()-1].P[2]; 
        // z of first sphere of first axon
        min_limits = first_axon.spheres[0].P[2]; 
        cout << "min_limits : " << min_limits << endl;
        cout << "max_limits : " << max_limits << endl;

        // voxel
        params.max_limits = Eigen::Vector3d(max_limits,max_limits,max_limits);
        params.min_limits = Eigen::Vector3d(min_limits,min_limits,min_limits);
        double volume = (params.max_limits[0]-params.min_limits[0])*(params.max_limits[1]-params.min_limits[1])*(params.max_limits[2]-params.min_limits[2]);
        //pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(min_limits,max_limits);
        //params.voxels_list.push_back(voxel_);
    
        cout << " Volume :" << volume << endl;
        // set icvf
        double icvf = computeICVF(params.min_limits, params.max_limits, dynamicsEngine->axons_list);
        // set number of particles
        if (params.ini_walker_flag == "intra"){
            params.setNumWalkers(params.concentration*volume*icvf);
        }
        else if (params.ini_walker_flag == "extra"){
            params.setNumWalkers(params.concentration*volume*(1.0-icvf));
        }
        else{
            params.setNumWalkers(params.concentration*volume);
        } 
        cout << "params.ini_walker_flag :" << params.ini_walker_flag << endl;
        
        cout << " ICVF :" << icvf<< endl;
        cout << " Number of particles :" << params.num_walkers << endl;
        cout << "voxel size :"  << max_limits << endl;



        in.close();
        
    }
}

void MCSimulation::addCylindersObstaclesFromFiles()
{

    for(unsigned i = 0; i < params.cylinders_files.size(); i++){


        std::ifstream in(params.cylinders_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}

            std::vector<std::string> jkr = split(line,' ');
            if (jkr.size() != 8){
                //std::cout << "\033[1;33m[Warning]\033[0m Cylinder orientation was set towards the Z direction by default" << std::endl;
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

        double x,y,z,r;
        double last_z = 0.0;
        int cyl_id, sph_id;
        int last_cyl_id = -1;
        std::string type_object;
        std::string header;
        int line_num = 0;

        Cylinder cylinder;

        int header_size = 8;

        for(unsigned j = 0; j < header_size; j++){  
            in >>header;
            //cout << "header :" << header << endl;
        } 

        while (in >>cyl_id >> sph_id >> type_object >> x >> y >> z >> r){

            // convert to mm
            x = x/1000.0;
            y = y/1000.0;
            z = z/1000.0;
            r = r/1000.0;
            if (line_num !=0 and last_cyl_id != cyl_id){
                // Calculate the factor to shift the decimal places
                double factor = std::pow(10, 2);
                // Round the number to the desired digits
                double rounded_last_z = std::round(last_z * factor) / factor;
                cylinder = Cylinder(cyl_id, Eigen::Vector3d(x,y,0), Eigen::Vector3d(x,y,rounded_last_z), r);
                    
                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  

                cylinder.setPercolation(perm_);

                // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                diff_i = params.diffusivity_intra; 
                diff_e = params.diffusivity_extra;
                cylinder.setDiffusion(diff_i, diff_e);

                dynamicsEngine->cylinders_list.push_back(cylinder);
            }
            last_z = z;
            last_cyl_id = cyl_id;
            line_num += 1;
        }
        in_perm.close();
        in.close();

        double max_limits = z; 
        // z of first sphere of first axon
        double min_limits = 0.0; 

        params.max_limits = Eigen::Vector3d(max_limits,max_limits,max_limits);
        params.min_limits = Eigen::Vector3d(min_limits,min_limits,min_limits);
        //pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(min_limits,max_limits);
        //params.voxels_list.push_back(voxel_);

        double volume = (params.max_limits[0]-params.min_limits[0])*(params.max_limits[1]-params.min_limits[1])*(params.max_limits[2]-params.min_limits[2]);
        
        std::cout << " Volume :" << volume << endl;
        // set icvf
        double icvf = computeAreaICVF(params.min_limits, params.max_limits, dynamicsEngine->cylinders_list);

        // set number of particles
        if (params.ini_walker_flag == "intra"){
            params.setNumWalkers(params.concentration*volume*icvf);
        }
        else if (params.ini_walker_flag == "extra"){
            params.setNumWalkers(params.concentration*volume*(1.0-icvf));
        }
        else{
            params.setNumWalkers(params.concentration*volume);
        } 
        std::cout << "params.ini_walker_flag :" << params.ini_walker_flag << endl;

        std::cout << " ICVF :" << icvf<< endl;
        std::cout << " Number of particles :" << params.num_walkers << endl;
        std::cout << "voxel size :"  << max_limits << endl;
        
    }
    /*
    for(unsigned i = 0; i < params.cylinders_files.size(); i++){

        bool z_flag = false;
        std::ifstream in(params.cylinders_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}

            std::vector<std::string> jkr = split(line,' ');
            if (jkr.size() != 7){
                z_flag = true;
                //std::cout << "\033[1;33m[Warning]\033[0m Cylinder orientation was set towards the Z direction by default" << std::endl;
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


        if(z_flag){
            double x,y,z,r;
            double scale;
            in >> scale;
            while (in >> x >> y >> z >> r)
            {
                Cylinder cyl(Cylinder(Eigen::Vector3d(x,y,z),Eigen::Vector3d(x,y,z+1.0),r,scale, perm_));

                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  
        
                cyl.setPercolation(perm_);

                // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                diff_i = params.diffusivity_intra; 
                diff_e = params.diffusivity_extra;
                cyl.setDiffusion(diff_i, diff_e);

                dynamicsEngine->cylinders_list.push_back(cyl);
            }
            in_perm.close();
            in.close();
        }
        else{
            double x,y,z,ox,oy,oz,r;
            double scale;
            in >> scale;
            while (in >> x >> y >> z >> ox >> oy >> oz >> r)
            {

                Cylinder cyl(Eigen::Vector3d(x,y,z),Eigen::Vector3d(ox,oy,oz),r,scale, perm_);

                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  
        
                cyl.setPercolation(perm_);

                // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                diff_i = params.diffusivity_intra; 
                diff_e = params.diffusivity_extra;
                cyl.setDiffusion(diff_i, diff_e);

                dynamicsEngine->cylinders_list.push_back(cyl);
            }
            in_perm.close();
            in.close();
        }

    }
    */
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


