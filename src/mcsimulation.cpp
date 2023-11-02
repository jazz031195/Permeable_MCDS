#include "mcsimulation.h"
#include <Eigen/Dense>
#include "simerrno.h"
#include "pgsesequence.h"
#include "gradientwaveform.h"
#include <iostream>
#include "constants.h"

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
                double l = (axons[i].spheres[j - 1].center - axons[i].spheres[j].center).norm(); // distance between centers
                double mean_r = (axons[i].spheres[j - 1].radius + axons[i].spheres[j].radius) / 2;

                if (withinBounds(min_limits, max_limits,axons[i].spheres[j].center, axons[i].spheres[j].radius) && withinBounds(min_limits, max_limits,axons[i].spheres[j-1].center, axons[i].spheres[j-1].radius))
                {
                    AreaC += l * M_PI * mean_r * mean_r;
                }
                else if (withinBounds(min_limits, max_limits,axons[i].spheres[j].center, 0) && withinBounds(min_limits, max_limits,axons[i].spheres[j-1].center, 0))
                {
                    AreaC += l * M_PI * mean_r * mean_r/2;
                }
            }
        }
    }
    return AreaC / AreaV; // ( total axons volume / total volume )
}

double computeICVF(Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector <Neuron> neurons)
{
    if (neurons.size() == 0)
        return 0;
    double AreaV = (max_limits[0] - min_limits[0]) * (max_limits[1] - min_limits[1]) * (max_limits[2] - min_limits[2]); // total volume
    double AreaC = 0;

    for (uint i = 0; i < neurons.size(); i++) // for all axons
    {
        for (uint d = 0; d < neurons[i].dendrites.size(); d++) // for all axons
        {
            for (uint s = 0; s < neurons[i].dendrites[d].subbranches.size(); s++) // for all axons
            {
                auto spheres = neurons[i].dendrites[d].subbranches[s].spheres;
                if (spheres.size() > 1)
                {
                    for (uint j = 1; j < spheres.size(); j++)
                    {
                        double l = (spheres[j - 1].center - spheres[j].center).norm(); // distance between centers
                        double mean_r = (spheres[j - 1].radius + spheres[j].radius) / 2;

                        if (withinBounds(min_limits, max_limits, spheres[j].center, spheres[j].radius) && withinBounds(min_limits, max_limits, spheres[j-1].center, spheres[j-1].radius))
                        {
                            AreaC += l * M_PI * mean_r * mean_r;
                        }
                        else if (withinBounds(min_limits, max_limits, spheres[j].center, 0) && withinBounds(min_limits, max_limits, spheres[j-1].center, 0))
                        {
                            AreaC += l * M_PI * mean_r * mean_r/2;
                        }
                    }
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

        double x,y,z,r,p;
        int ax_id, sph_id;
        int last_ax_id = -1;
        std::string type_object;
        std::string header;

        std::vector<Sphere> spheres_ ;
        Sphere sphere_;
        int line_num = 0;

        int header_size = 8;

        for(int j = 0; j < header_size; j++)
            in >> header;
        

        while (in >> ax_id >> sph_id >> type_object >> x >> y >> z >> r >> p){

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
        max_limits = first_axon.spheres[first_axon.spheres.size()-1].center[2]; 
        // z of first sphere of first axon
        min_limits = first_axon.spheres[0].center[2]; 
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
        if (params.concentration != 0){
            if (params.ini_walker_flag == "intra"){
                params.setNumWalkers(params.concentration*volume*icvf);
            }
            else if (params.ini_walker_flag == "extra"){
                params.setNumWalkers(params.concentration*volume*(1.0-icvf));
            }
            else{
                params.setNumWalkers(params.concentration*volume);
            } 
        }
        cout << "params.ini_walker_flag :" << params.ini_walker_flag << endl;
        
        cout << " ICVF :" << icvf<< endl;
        cout << " Number of particles :" << params.num_walkers << endl;
        cout << "voxel size :"  << max_limits << endl;



        in.close();
        
    }
}

void MCSimulation::printSubstrate(ostream &out) 
{
    out << 1 << endl; //scale
    out << 0 << endl; //volume_inc_perc
    out << 0 << endl; //dyn_perc
    out << dynamicsEngine->icvf << endl;
    out << dynamicsEngine->params.min_limits[0] << endl; //min_limits [mm]
    out << dynamicsEngine->params.max_limits[0] << endl; //max_limits [mm]

    auto neurons = dynamicsEngine->neurons_list;
    for (unsigned i = 0; i < neurons.size(); i++)
    {
        // Print for soma : x y z r bool_active
        out << neurons[i].soma.center[0] << " " 
        << neurons[i].soma.center[1] << " "
        << neurons[i].soma.center[2] << " "
        << neurons[i].soma.radius    << endl; 
        out << "Soma " + to_string(i) << endl;

        for (size_t j = 0; j < neurons[i].dendrites.size(); j++)
        {
            for (size_t k = 0; k < neurons[i].dendrites[j].subbranches.size(); k++)
            {
                for (size_t l = 0; l < neurons[i].dendrites[j].subbranches[k].spheres.size(); l++)
                {
                    // Print for each dendrite, each sphere
                    out << neurons[i].dendrites[j].subbranches[k].spheres[l].center[0] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[1] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[2] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].radius << endl; 
                        
                }
                out << "Segment " + to_string(k);
                if(neurons[i].dendrites[j].subbranches[k].proximal_branching.size() == 2)
                    out << " proximal " + to_string(neurons[i].dendrites[j].subbranches[k].proximal_branching[0]) + " " + to_string(neurons[i].dendrites[j].subbranches[k].proximal_branching[1]);
                else
                    out << " proximal " + to_string(neurons[i].dendrites[j].subbranches[k].proximal_branching[0]);

                if(neurons[i].dendrites[j].subbranches[k].distal_branching.size() == 2)
                    out << " distal " + to_string(neurons[i].dendrites[j].subbranches[k].distal_branching[0]) + " " + to_string(neurons[i].dendrites[j].subbranches[k].distal_branching[1]) << endl;
                else
                    out << " distal " + to_string(neurons[i].dendrites[j].subbranches[k].distal_branching[0]) << endl;

            }
            out << "Dendrite " + to_string(j) << endl;
        }
        out << "Neuron " + to_string(i) << endl;
    }
}

void MCSimulation::readNeurons_fromSWC(int const& neurons_files_id)
{
    std::ifstream in(params.neurons_files[neurons_files_id]);

    if(!in){
        std::cout <<  "[ERROR] Unable to open:" << params.neurons_files[neurons_files_id] << std::endl;
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
        if (jkr.size() != 7){
            std::cout << jkr.size() <<  " elements per line" << std::endl;
            std::cout << "wrong number of elements per line in file" << std::endl;
        }
        break;
    }
    in.close();

    // Permeability file - if any
    double perm_; 

    std::ifstream in_perm;
    if(params.neuron_permeability_files.size() > 0)
        in_perm.open(params.neuron_permeability_files[neurons_files_id]);
    

    // Local permeability - Different for each obstacle
    if(params.neuron_permeability_files.size() > 0)
        in_perm >> perm_;
    // Global permeability - Same for all obstacle
    else
        perm_ = params.obstacle_permeability;
    
            
    // Diffusion coefficients
    double diff_i = params.diffusivity_intra; 
    double diff_e = params.diffusivity_extra;

    in.open(params.neurons_files[neurons_files_id]);
    // cout << params.neurons_files[neurons_files_id] << endl;
    double l, x, y, z, r, p;
    double scale;
    double volume_inc_perc, dyn_perc, icvf;
    double max_limits, min_limits;

    in >> scale;
    in >> volume_inc_perc;
    in >> dyn_perc;
    in >> icvf;
    in >> min_limits;
    in >> max_limits;

    std::vector<Sphere> spheres_ ;
    std::vector<Dendrite> dendrites_ ;
    Sphere sphere_;
    Sphere soma;
    string part;
    int id = 0;
    vector<vector<double>> lines;
    std::string line; getline( in, line );
    for( std::string line; getline( in, line ); ){
        std::vector<std::string> jkr = split(line,' ');
        int ax_id = 0;
        l = stod(jkr[0]);
        x = stod(jkr[2]);
        y = stod(jkr[3]);
        z = stod(jkr[4]);
        r = stod(jkr[5]);
        p = stod(jkr[6]);

        // Soma
        if(p == -1)
        {
            soma = Sphere(id, ax_id, Eigen::Vector3d(x,y,z), r, scale);
            soma.setPercolation(perm_);
            soma.setDiffusion(diff_i, diff_e);
        }
        // Start dendrite
        else if(p == 1)
        {
            // End of dendrite => create it
            if(lines.size() > 0)
            {
                Dendrite dendrite_ = createDendrites(lines, scale, perm_);
                dendrite_.subbranches[0].spheres[0].add_neighbor(new Sphere(soma));
                soma.add_neighbor(new Sphere(dendrite_.subbranches[0].spheres[0]));
                dendrites_.push_back(dendrite_);
                
                // Start fresh for the next dendrite
                lines.clear();
            }
            lines.push_back({l, x, y, z, r, p});
        }
        else
            lines.push_back({l, x, y, z, r, p});
    } 

    // Add the last dendrite
    Dendrite dendrite_ = createDendrites(lines, scale, perm_);
    dendrite_.subbranches[0].spheres[0].add_neighbor(new Sphere(soma));
    soma.add_neighbor(new Sphere(dendrite_.subbranches[0].spheres[0]));
    dendrites_.push_back(dendrite_);
    lines.clear();

    Neuron neuron(dendrites_, soma);

    // Add the permeabilities to all the objects
    for(size_t d=0; d < neuron.dendrites.size(); ++d)
    {
        for(size_t s=0; s < neuron.dendrites[d].subbranches.size(); ++s)
        {
            for(size_t sph=0; sph < neuron.dendrites[d].subbranches[s].spheres.size(); ++sph)
            {
                neuron.dendrites[d].subbranches[s].spheres[sph].setPercolation(perm_);
                neuron.dendrites[d].subbranches[s].spheres[sph].setDiffusion(diff_i, diff_e);
            }
        }
    }
    neuron.add_projection();
    neuron.setPercolation(perm_);
    neuron.setDiffusion(diff_i, diff_e);
    dynamicsEngine->neurons_list.push_back(neuron);
    dynamicsEngine->area = dynamicsEngine->area + neuron.get_Area();

    // string file = params.output_base_name + "_neurons_list.txt";
    // ofstream out(file);
    // printSubstrate(out); 

    dendrites_.clear();       
}



/* Function that finds the proximal indices to a segment 
   \param segments (vector<vector<double>>) : all the segments of a dendrite
   \param segment_id (int) : indice of the segment of interest
*/
vector<int> find_proximal(vector<vector<double>> const& segments, int const& segment_id)
{x
    vector<int> ids;
    // If first segment => linked to soma => id proximal = -1
    if(segments[segment_id][0] == 0)
        return {-1};

    for (int i=0; i < segments.size(); ++i) 
    {
        if((segments[i][1] == segments[segment_id][0]) && (i != segment_id))
            ids.push_back(i);
        if((segments[i][0] == segments[segment_id][0]) && (i != segment_id))
            ids.push_back(i);
    }
    return ids; 
}

/* Function that finds the distal indices to a segment 
   \param segments (vector<vector<double>>) : all the segments of a dendrite
   \param segment_id (int) : indice of the segment of interest
*/
vector<int> find_distal(std::vector<vector<double>> const& segments, int const& segment_id)
{
    vector<int> ids;

    for (int i=0; i < segments.size(); ++i) 
    {
        if(segments[i][0] == segments[segment_id][1])
            ids.push_back(i);
    }

    if(ids.size() > 0)
        return ids; 
    else
        return {-1};
}

Dendrite MCSimulation::createDendrites(vector<vector<double>> const& lines, double const& scale, double& perm_)
{
    std::vector<Axon> subbranches_ ;
    Dendrite dendrite;

    std::vector<vector<double>> segments;
    for(size_t i=0; i < lines.size(); ++i)
    {
        // When last entry == 1, it is linked to the soma
        if(lines[i][lines[i].size() - 1] == 1)
            continue;

        // Shift by which to shift the parent id so that it starts at 0
        int shift_id  = lines[1][lines[1].size() - 1];  
        int parent_id = lines[i][lines[i].size() - 1] - shift_id;  
        segments.push_back({lines[parent_id][0] - shift_id, lines[i][0] - shift_id});
    }

    // Sort the 2D vector based on the first column
    std::sort(segments.begin(), segments.end(), [](const std::vector<double>& a, const std::vector<double>& b) 
    {
        if (a[0] != b[0])
            return a[0] < b[0];
        // If the 1st column is at equality, check the second
        else
            return a[1] < b[1];
    });

    Eigen::Vector3d parent, child, begin;
    double radius, funnel_radius, target_radius;
    Eigen::Vector3d dir;
    vector<Sphere> spheres_to_add;
    // Display the sorted 2D vector
    for (size_t i=0; i < segments.size(); ++i) {

        parent = {lines[segments[i][0]][1], lines[segments[i][0]][2], lines[segments[i][0]][3]};
        child  = {lines[segments[i][1]][1], lines[segments[i][1]][2], lines[segments[i][1]][3]};
        radius = target_radius = lines[segments[i][1]][4];
        dir    = (child - parent).normalized();
        
        vector<int> id_prox = find_proximal(segments, i);
        vector<int> id_dist = find_distal(segments, i);

        Eigen::Vector3d center = parent;
        int j = 0;
        double funnel_radius = radius;
        if(params.funnel)
            funnel_radius = 1.5e-3; // TODO [ines] : find optimal

        double slope = 1.0 / 5.0; // radius from 1.5 to 0.5 over 5um
        // Continue until the next center reached the child center
        while((center - child).norm() > 1e-5)
        {            
            // If 1st segment => starts from 0
            if(i == 0)
            {
                radius = funnel_radius - slope * double(j) * (radius / params.sphere_overlap);
                if(radius <= target_radius)
                    radius = target_radius;
                center = parent + double(j) * (target_radius / params.sphere_overlap) * dir;
            }
            // If not => starts from 1, so that no overlay of spheres at the branching
            else
                center = parent + double(j+1) * (radius / params.sphere_overlap) * dir;
    
            Sphere s(j, i, center, radius, scale);
            if(j > 0)
            {
                s.add_neighbor(new Sphere(spheres_to_add[spheres_to_add.size() - 1]));
                spheres_to_add[spheres_to_add.size() - 1].add_neighbor(new Sphere(s));
            }
            spheres_to_add.push_back(s);
            j++;
        }

        Axon subbranch(i, radius, begin, begin, id_prox, id_dist);
        subbranch.set_spheres(spheres_to_add);
        dendrite.add_subbranch(subbranch);
        spheres_to_add.clear();
    }

    for(size_t j=0; j < dendrite.subbranches.size(); j++)
    {
        vector<int> dist = dendrite.subbranches[j].distal_branching;
        if(dist.size() == 2)
        {
            dendrite.subbranches[j].spheres[dendrite.subbranches[j].spheres.size() - 1].add_neighbor(new Sphere(dendrite.subbranches[dist[0]].spheres[0]));
            dendrite.subbranches[dist[0]].spheres[0].add_neighbor(new Sphere(dendrite.subbranches[j].spheres[dendrite.subbranches[j].spheres.size() - 1]));
            dendrite.subbranches[dist[0]].spheres[0].add_neighbor(new Sphere(dendrite.subbranches[dist[1]].spheres[0]));
            dendrite.subbranches[j].spheres[dendrite.subbranches[j].spheres.size() - 1].add_neighbor(new Sphere(dendrite.subbranches[dist[1]].spheres[0]));
            dendrite.subbranches[dist[1]].spheres[0].add_neighbor(new Sphere(dendrite.subbranches[j].spheres[dendrite.subbranches[j].spheres.size() - 1]));
            dendrite.subbranches[dist[1]].spheres[0].add_neighbor(new Sphere(dendrite.subbranches[dist[0]].spheres[0]));    
        }
    }
    dendrite.add_projection();
    return dendrite;
}

void MCSimulation::readNeurons_fromList(int const& neurons_files_id)
{
    std::ifstream in(params.neurons_files[neurons_files_id]);

    if(!in){
        std::cout <<  "[ERROR] Unable to open:" << params.neurons_files[neurons_files_id] << std::endl;
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
    if(params.neuron_permeability_files.size() > 0)
        in_perm.open(params.neuron_permeability_files[neurons_files_id]);
    

    // Local permeability - Different for each obstacle
    if(params.neuron_permeability_files.size() > 0)
        in_perm >> perm_;
    // Global permeability - Same for all obstacle
    else
        perm_ = params.obstacle_permeability;
    
            
    // Diffusion coefficients
    double diff_i; 
    double diff_e;

    in.open(params.neurons_files[neurons_files_id]);
    // cout << params.neurons_files[neurons_files_id] << endl;
    double x,y,z,r;
    double scale;
    double volume_inc_perc, dyn_perc, icvf;
    double max_limits, min_limits;

    in >> scale;
    in >> volume_inc_perc;
    in >> dyn_perc;
    in >> icvf;
    in >> min_limits;
    in >> max_limits;

    std::vector<Sphere> spheres_ ;
    std::vector<Axon> subbranches_ ;
    std::vector<Dendrite> dendrites_ ;
    Sphere sphere_;
    Sphere soma;
    string part;
    int id;

    for( std::string line; getline( in, line ); ){
                
        std::vector<std::string> jkr = split(line,' ');
        int ax_id = 0;
        if(jkr.size() == 4){
            x = stod(jkr[0]);
            y = stod(jkr[1]);
            z = stod(jkr[2]);
            r = stod(jkr[3]);
            sphere_ = Sphere(id, ax_id, Eigen::Vector3d(x,y,z), r, scale);
            if(spheres_.size() > 0)
            {
                sphere_.add_neighbor(new Sphere(spheres_[spheres_.size() - 1]));
                spheres_[spheres_.size() - 1].add_neighbor(new Sphere(sphere_));
            }
            
            spheres_.push_back(sphere_);
    
            if(subbranches_.size() > 0 && spheres_.size() == 0)
            {
                sphere_.add_neighbor(new Sphere(subbranches_[0].spheres[0]));
                subbranches_[0].spheres[0].add_neighbor(new Sphere(sphere_));
                
            }
            else if (spheres_.size() == 0)
            {
                soma.add_neighbor(new Sphere(sphere_));
                sphere_.add_neighbor(new Sphere(soma));
            }
            // cout << "adding sphere, radius :" << sphere_.radius  << endl;
        
        }
        if(jkr.size() == 2){
            part = jkr[0];
            id = stod(jkr[1]);

            // If flag "soma", create a soma
            if( part.find("Soma") != std::string::npos && spheres_.size() == 1)
            {
                soma = spheres_[0];
                soma.setPercolation(perm_);
                // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                diff_i = params.diffusivity_intra; 
                diff_e = params.diffusivity_extra;
                soma.setDiffusion(diff_i, diff_e);
                spheres_.clear();
            }
            // If dendrite, create it from spheres_ and store it into axons_
            else if( part.find("Dendrite") != std::string::npos && subbranches_.size() > 0)
            {
                Dendrite dendrite;

                for(size_t i=0; i < subbranches_.size(); ++i)
                {
                    vector<int> prox_id = subbranches_[i].proximal_branching;
                    if(prox_id.size() == 2)
                    {
                        subbranches_[i].spheres[0].add_neighbor(new Sphere(subbranches_[prox_id[0]].spheres[subbranches_[prox_id[0]].spheres.size() - 1]));
                        subbranches_[i].spheres[0].add_neighbor(new Sphere(subbranches_[prox_id[1]].spheres[subbranches_[prox_id[1]].spheres.size() - 1]));
                    }

                    vector<int> dist_id = subbranches_[i].distal_branching;
                    if(dist_id.size() == 2)
                    {
                        subbranches_[i].spheres[subbranches_[i].spheres.size() - 1].add_neighbor(new Sphere(subbranches_[dist_id[0]].spheres[0]));
                        subbranches_[i].spheres[subbranches_[i].spheres.size() - 1].add_neighbor(new Sphere(subbranches_[dist_id[1]].spheres[0]));
                    }
                }

                dendrite.set_dendrite(subbranches_);               
                dendrites_.push_back(dendrite);
                subbranches_.clear();
                // cout << "adding dendrite : "  << id << endl;
            }
            // If neuron, create it from soma and axons
            else if( part.find("Neuron") != std::string::npos)
            {
                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                } 
                Neuron neuron(dendrites_, soma);

                for(size_t d=0; d < neuron.dendrites.size(); ++d)
                {
                    for(size_t s=0; s < neuron.dendrites[d].subbranches.size(); ++s)
                    {
                        for(size_t sph=0; sph < neuron.dendrites[d].subbranches[s].spheres.size(); ++sph)
                        {
                            neuron.dendrites[d].subbranches[s].spheres[sph].setPercolation(perm_);
                            // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                            diff_i = params.diffusivity_intra; 
                            diff_e = params.diffusivity_extra;
                            neuron.dendrites[d].subbranches[s].spheres[sph].setDiffusion(diff_i, diff_e);
                        }
                    }
                }
                neuron.add_projection();
                neuron.setPercolation(perm_);
                neuron.setDiffusion(diff_i, diff_e);
                dynamicsEngine->neurons_list.push_back(neuron);
                dynamicsEngine->area = dynamicsEngine->area + neuron.get_Area();
                dendrites_.clear();
                // TODO : why nb_dendrites not printed ? [ines]
                // cout << "adding neuron: "  << id << ", nb_dendrites: " << neuron.nb_dendrites << endl;
            }
        }
        if(jkr.size() > 4)
        {
            part = jkr[0];
            id   = stod(jkr[1]);
            auto it = std::find_if( std::begin( jkr ), std::end( jkr ),
                            [&]( const string s ){ return s == "distal"; } );
        
            const int pos = std::distance( jkr.begin(), it );

            vector<int> prox_id = {int(stod(jkr[3]))};
            if((pos - 3) == 2)
                prox_id.push_back(int(stod(jkr[4])));
            
            vector<int> dist_id = {int(stod(jkr[pos + 1]))};
            if((pos + 1) != (jkr.size() - 1)) // TODO [ines] : check if parentheses needed
                dist_id.push_back(int(stod(jkr[pos + 2])));

            // If Segment, create it from spheres_ and store it into axons_
            if( part.find("Segment") != std::string::npos && spheres_.size() > 0)
            {
                Eigen::Vector3d begin = {min_limits, min_limits, min_limits};
                Eigen::Vector3d end   = {max_limits, max_limits, max_limits};
                Axon subbranch(id, r, begin, end, prox_id, dist_id);
                subbranch.set_spheres(spheres_);

                for (unsigned i = 0; i < spheres_.size(); i++){
                    spheres_[i].setPercolation(perm_);
                    // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                    diff_i = params.diffusivity_intra; 
                    diff_e = params.diffusivity_extra;
                    spheres_[i].setDiffusion(diff_i, diff_e);
                }

                spheres_.clear();
                subbranches_.push_back(subbranch);
                // cout << "adding segment "  << endl;
                //TODO [ines] : add proximal & distal branching reading
            }
        }
    }

    params.max_limits = Eigen::Vector3d(max_limits, max_limits, max_limits);
    params.min_limits = Eigen::Vector3d(min_limits, min_limits, min_limits);
    params.gamma_icvf = icvf;

    in.close();

    double volume          = (params.max_limits[0] - params.min_limits[0]) * (params.max_limits[1] - params.min_limits[1]) * (params.max_limits[2] - params.min_limits[2]);
    // double icvf_calculated = computeICVF(params.min_limits, params.max_limits, dynamicsEngine->neurons_list);

    if (params.concentration != 0){
        if (params.ini_walker_flag == "intra")
            params.setNumWalkers(params.concentration * volume * icvf);
        else if (params.ini_walker_flag == "extra")
            params.setNumWalkers(params.concentration * volume * (1.0 - icvf));
        else
            params.setNumWalkers(params.concentration * volume);
        
    }
}

void MCSimulation::addNeuronsObstaclesFromFiles()
{

    // Read neurons from file
    for(unsigned i = 0; i < params.neurons_files.size(); i++){

        if(params.neurons_files[i].find(".swc") != std::string::npos)
            readNeurons_fromSWC(i);
        else
            readNeurons_fromList(i);
        
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

        for(int j = 0; j < header_size; j++) 
            in >> header;
        

        while (in >> cyl_id >> sph_id >> type_object >> x >> y >> z >> r){

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
        if (params.concentration != 0){
            if (params.ini_walker_flag == "intra"){
                params.setNumWalkers(params.concentration*volume*icvf);
            }
            else if (params.ini_walker_flag == "extra"){
                params.setNumWalkers(params.concentration*volume*(1.0-icvf));
            }
            else{
                params.setNumWalkers(params.concentration*volume);
            } 
        }
        std::cout << "params.ini_walker_flag :" << params.ini_walker_flag << endl;

        std::cout << " ICVF :" << icvf << endl;
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


