#include "parallelmcsimulation.h"
#include <iomanip>
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include "gammadistribution.h"
#include "uniformdistribution.h"
#include "gaussiandistribution.h"
#include "spheredistribution.h"
#include "cylinderdistribution.h"


using namespace std;

ParallelMCSimulation::ParallelMCSimulation(std::string config_file)
{
    mean_second_passed = 0;

    total_sim_particles = 0;

    SimErrno::checkConfigurationFile(config_file.c_str());

    cout << "read scheme files" << endl;

    params.readSchemeFile(config_file);

    cout << "read scheme files" << endl;

    //printSimulationInfo();

    SimErrno::checkSimulationParameters(params);
    SimErrno::printSimulatinInfo(params,std::cout);

    initializeUnitSimulations();
}

ParallelMCSimulation::ParallelMCSimulation(Parameters &params)
{
    this->params = params;
    mean_second_passed = 0;
    total_sim_particles = 0;
    SimErrno::checkSimulationParameters(params);
    SimErrno::printSimulatinInfo(params,std::cout);
    this->params = params;
    initializeUnitSimulations();
    icvf=0;
}

ParallelMCSimulation::~ParallelMCSimulation()
{
    for(unsigned i = 0; i < simulations.size(); i++){
        delete simulations[i];
    }
}

void ParallelMCSimulation::startSimulation()
{
    cout<<setfill('-');
    cout << SH_FG_PURPLE << "/********************   MC/DC Simulation START:  *************************/" << SH_DEFAULT << "\n";

    unsigned int count_simu=0, max_nb_proc=50;

    for(unsigned int i =0; i < simulations.size(); i++){
    
        sim_threads.push_back(std::thread(&MCSimulation::startSimulation,(simulations[i])));
        count_simu++;

        if((count_simu==max_nb_proc) || (i==simulations.size()-1)){
            for(unsigned j=i-count_simu+1; j <= i; j++){
                sim_threads[j].join();
            }
            count_simu = 0;
        }
    }

    /*
    for(unsigned i =0; i < simulations.size(); i++){
        sim_threads[i].join();
    }
    */
    for(unsigned i =0; i < simulations.size(); i++){
        mean_second_passed  += simulations[i]->dynamicsEngine->second_passed/double(simulations.size());
        total_sim_particles += simulations[i]->dynamicsEngine->num_simulated_walkers;
    }

    cout<<setfill('-');

    SimErrno::info( " All " + to_string(params.num_proc) +  " simulations ended after: " + DynamicsSimulation::secondsToMinutes(mean_second_passed)
                    + " in average",cout  );

    SimErrno::info( " Joining resulting data... ",cout  );

    jointResults();

    SimErrno::info( " Done.",cout  );

    ofstream out(params.output_base_name+"_simulation_info.txt", std::ofstream::app);
    SimErrno::info("All " + to_string(params.num_proc) +  " simulations ended after: " + DynamicsSimulation::secondsToMinutes(mean_second_passed)
                   + " in average",out,false);
    SimErrno::info("Number of particles labeled as stuck: "        + to_string(stuck_count)  ,out,false);
    SimErrno::info("Number of particles eliminated due crossings: "+ to_string(illegal_count),out,false);
    int tot_nbr_walker_extra_ini = 0;
    int tot_nbr_walker_intra_ini  = 0;
    int tot_nbr_walker_axons_ini  = 0;
    int tot_nbr_walker_glials_ini  = 0;
    int tot_nbr_walker_extra_final = 0;
    int tot_nbr_walker_intra_final  = 0;
    int tot_nbr_walker_axons_final  = 0;
    int tot_nbr_walker_glials_final  = 0;

    int tot_nbr_bounces = 0;
    int tot_nbr_legal_crossings = 0;
    for (unsigned i = 0; i < simulations.size(); i++){
        tot_nbr_walker_extra_ini += simulations[i]->dynamicsEngine->nbr_walker_extra_ini;
        tot_nbr_walker_intra_ini += simulations[i]->dynamicsEngine->nbr_walker_intra_ini;
        tot_nbr_walker_axons_ini += simulations[i]->dynamicsEngine->nbr_walker_axons_ini;
        tot_nbr_walker_glials_ini += simulations[i]->dynamicsEngine->nbr_walker_glials_ini;
        tot_nbr_walker_extra_final += simulations[i]->dynamicsEngine->nbr_walker_extra_final;
        tot_nbr_walker_intra_final += simulations[i]->dynamicsEngine->nbr_walker_intra_final;
        tot_nbr_walker_axons_final += simulations[i]->dynamicsEngine->nbr_walker_axons_final;
        tot_nbr_walker_glials_final += simulations[i]->dynamicsEngine->nbr_walker_glials_final;

        tot_nbr_bounces += simulations[i]->dynamicsEngine->tot_nbr_bounces;
        tot_nbr_legal_crossings += simulations[i]->dynamicsEngine->tot_nbr_legal_crossings;

    }
    SimErrno::info("Number of intracellular particles at the start : "+ to_string(tot_nbr_walker_intra_ini),out,false);
    SimErrno::info("Number of extracellular particles at the start : "+ to_string(tot_nbr_walker_extra_ini),out,false);
    SimErrno::info("Number of axonal particles at the start : "+ to_string(tot_nbr_walker_axons_ini),out,false);
    SimErrno::info("Number of glial particles at the start : "+ to_string(tot_nbr_walker_glials_ini),out,false);
    SimErrno::info("Number of intracellular particles at the end : "+ to_string(tot_nbr_walker_intra_final),out,false);
    SimErrno::info("Number of extracellular particles at the end : "+ to_string(tot_nbr_walker_extra_final),out,false);
    SimErrno::info("Number of walkers in axons at the end : "+ to_string(tot_nbr_walker_axons_final),out,false);
    SimErrno::info("Number of walkers in glial cells at the end : "+ to_string(tot_nbr_walker_glials_final),out,false);

    SimErrno::info("Total number of bounces: "                     + to_string(tot_nbr_bounces),out,false);
    SimErrno::info("Total number of legal crossings: "             + to_string(tot_nbr_legal_crossings),out,false);
    
    if(params.max_simulation_time > 0){
        SimErrno::info("Number of simulated particles: "+ to_string(total_sim_particles),out,false);
        SimErrno::info("Mean simulation speed: "+ to_string( unsigned(params.num_steps*total_sim_particles/mean_second_passed)) + " steps/second",out,false);
    }
    else{
        SimErrno::info("Mean simulation speed: "+ to_string( unsigned(params.num_steps*params.num_walkers/mean_second_passed)) + " steps/second",out,false);
    }
    if(params.voxels_list.size()>0){
        string message = "Voxel limits:";
        SimErrno::info(message,out,false);

        out << std::setprecision(10) << "( " <<  params.voxels_list[0].first[0] <<
                " " << params.voxels_list[0].first[1]<<
                " " << params.voxels_list[0].first[2]<< " )\n";

        out << std::setprecision(10) << "( " <<  params.voxels_list[0].second[0]<<
                " " << params.voxels_list[0].second[1] <<
                " " << params.voxels_list[0].second[2] << " )\n";
    }

    if(params.custom_sampling_area){
        out << std::fixed;
        string message="Custom sampled area (mm):";
        SimErrno::info(message,out,false);
        out << std::setprecision(10) << "( " <<  params.min_sampling_area[0] <<
                " " << params.min_sampling_area[1]<<
                " " << params.min_sampling_area[2]<< " )\n";

        out << std::setprecision(10) << "( " <<  params.max_sampling_area[0]<<
                " " << params.max_sampling_area[1] <<
                " " << params.max_sampling_area[2] << " )\n";
    }

    if((params.custom_sampling_area or params.voxels_list.size()>0) and (params.computeVolume)){
        string message="Estimated Intra-axonal volume from sampling (mm^2):";
        out << std::scientific;
        SimErrno::info(message,out,false);
        out << std::scientific;
        out << aprox_volumen << endl;
    }
    out.close();
}


void ParallelMCSimulation::initializeUnitSimulations()
{

    //Build anything that needs to be syn between simulations.
    specialInitializations();
   

    // The number of walker N devided between all the processes
    unsigned N_per_sim = params.num_walkers/params.num_proc;

    for(unsigned i = 0; i < params.num_proc-1; i++){

        //Parameters for each simulation simulation
        Parameters params_temp       = params;

        params_temp.num_walkers      = N_per_sim;

        params_temp.output_base_name+= "_"+std::to_string(i);

        MCSimulation* simulation_ = new MCSimulation(params_temp);

        
        simulation_->dynamicsEngine->print_expected_time = 0;
        simulations.push_back(simulation_);

        if(params.verbatim)
            SimErrno::info( " Sim: " + to_string(simulation_->dynamicsEngine->id) + " Initialized",cout);
    }

    //We set thet remaining walkers so it sums up the desired total

    //Parameters for each simulation simulation
    Parameters params_temp = params;
    params_temp.num_walkers = params.num_walkers - N_per_sim *(params.num_proc-1);
    params_temp.output_base_name+= "_"+std::to_string(params.num_proc-1);

    MCSimulation* simulation_ = new MCSimulation(params_temp);
    simulations.push_back(simulation_);

    if(params.verbatim)
        SimErrno::info( " Sim: " + to_string(simulation_->dynamicsEngine->id) + " Initialized",cout);


}

void ParallelMCSimulation::jointResults()
{

    //====================================================================================================
    //Complete DWI data.
    //====================================================================================================

    //Writes the ouput data
    if(simulations[0]->dataSynth){

        if (simulations[0]->dataSynth->type != "PGSE_INTERVALS"){

            std::string outDWI   = params.output_base_name  + "_DWI.txt";
            std::string outDWIi   = params.output_base_name  + "_DWI_img.txt";
            std::string outPhase  = params.output_base_name  + "_phase_shift.txt";

            std::string boutDWI    = params.output_base_name  + "_DWI.bfloat";
            std::string boutDWIi   = params.output_base_name  + "_DWI_img.bfloat";
            std::string boutPhase  = params.output_base_name  + "_phase_shift.bfloat";

            std::ofstream phase_out,phase_bout, dwi_out, dwii_out, dwi_bout, dwii_bout;

            if(params.log_phase_shift){
                if(params.write_bin)
                    phase_bout.open(boutPhase, std::ofstream::binary);   //phase shift (binary)

                if(params.write_txt)
                    phase_out.open(outPhase, std::ofstream::out);     //phase shift  (txt)
            }

            if(params.write_txt){
                dwi_out.open(outDWI  ,std::ofstream::out);       //real part
                dwii_out.open(outDWIi,std::ofstream::out);       //img part
            }

            if(params.write_bin){
                dwi_bout.open(boutDWI  ,std::ofstream::binary);       //real part (binary)
                dwii_bout.open(boutDWIi,std::ofstream::binary);       //img part  (binary)
            }
            for (unsigned i = 0 ; i < simulations[0]->dataSynth->DWI.size(); i++ )
            {
                double DWI   = 0;
                double DWIi  = 0;

                std::vector<int> phase(3600, 0);
                for(unsigned p = 0; p < params.num_proc; p++)
                {
                    DWI+=  simulations[p]->dataSynth->DWI[i];       //real part
                    DWIi+= simulations[p]->dataSynth->DWIi[i];      //img  part

                    // for each histogram bin, 36000 fixed size
                    for(unsigned h = 0; h < 3600; h++)
                    {
                        phase[h]+=simulations[p]->dataSynth->phase_shift_distribution(i,h);
                    }

                }


                if(params.write_txt){
                    dwi_out  << DWI << std::endl;
                    dwii_out << DWIi << std::endl;
                }

                if(params.write_bin){
                    float holder = float(DWI);
                    dwi_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                    holder = float(DWIi);
                    dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                }

                if(params.log_phase_shift){

                    if(params.write_txt){
                        //write the full histogram for the gradient i
                        for(unsigned h = 0; h < 3600; h++)
                        {
                            phase_out << phase[h];
                            if(h < 3600-1)
                                phase_out << " ";
                            else
                                phase_out << endl;
                        }
                    }

                    if(params.write_bin){
                        //write the full histogram for the gradient i
                        for(unsigned h = 0; h < 3600; h++)
                        {
                            float holder = phase[h];
                            phase_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                        }
                    }
                }
            }
            if(params.write_txt){
                dwi_out.close();
                dwii_out.close();
            }
            if(params.write_bin){
                dwi_bout.close();
                dwii_bout.close();
            }

            if(params.log_phase_shift)
                phase_out.close();
        }
        else{

            std::string outDWI_intra   = params.output_base_name  + "_DWI_intra.txt";
            std::string outDWIi_intra   = params.output_base_name  + "_DWI_intra.txt";
            std::string outPhase_intra  = params.output_base_name  + "_phase_shift_intra.txt";

            std::string outDWI_extra   = params.output_base_name  + "_DWI_extra.txt";
            std::string outDWIi_extra   = params.output_base_name  + "_DWI_extra.txt";
            std::string outPhase_extra  = params.output_base_name  + "_phase_shift_extra.txt";

            std::string boutDWI_intra    = params.output_base_name  + "_DWI_intra.bfloat";
            std::string boutDWIi_intra   = params.output_base_name  + "_DWI_img_intra.bfloat";
            std::string boutPhase_intra  = params.output_base_name  + "_phase_shift_intra.bfloat";

            std::string boutDWI_extra    = params.output_base_name  + "_DWI_extra.bfloat";
            std::string boutDWIi_extra   = params.output_base_name  + "_DWI_img_extra.bfloat";
            std::string boutPhase_extra  = params.output_base_name  + "_phase_shift_extra.bfloat";

            std::string bout_compartment_walkers  = params.output_base_name  + "_walkers_compartment.bfloat";
            std::string out_compartment_walkers  = params.output_base_name  + "_walkers_compartment.txt";

            std::ofstream phase_out_intra, phase_out_extra, dwi_out_intra, dwii_out_extra, compartment_out, dwii_out_intra, dwi_out_extra;
            std::ofstream phase_bout_intra, phase_bout_extra,  dwi_bout_intra, dwii_bout_extra, compartment_bout, dwii_bout_intra, dwi_bout_extra;

            if(params.log_phase_shift){
                if(params.write_bin)
                    phase_bout_intra.open(boutPhase_intra, std::ofstream::binary);   //phase shift (binary)
                    phase_bout_extra.open(boutPhase_extra, std::ofstream::binary);   //phase shift (binary)

                if(params.write_txt)
                    phase_out_intra.open(outPhase_intra, std::ofstream::out);     //phase shift  (txt)
                    phase_out_extra.open(outPhase_extra, std::ofstream::out);     //phase shift  (txt)
            }

            if(params.write_txt){
                dwi_out_intra.open(outDWI_intra  ,std::ofstream::out);       //real part
                dwii_out_intra.open(outDWIi_intra,std::ofstream::out);       //img part

                dwi_out_extra.open(outDWI_extra  ,std::ofstream::out);       //real part
                dwii_out_extra.open(outDWIi_extra,std::ofstream::out);       //img part

                compartment_out.open(out_compartment_walkers, std::ofstream::out); //compartment walkers
            }

            if(params.write_bin){
                dwi_bout_intra.open(boutDWI_intra  ,std::ofstream::binary);       //real part (binary)
                dwii_bout_intra.open(boutDWIi_intra,std::ofstream::binary);       //img part  (binary)

                dwi_bout_extra.open(boutDWI_extra  ,std::ofstream::binary);       //real part (binary)
                dwii_bout_extra.open(boutDWIi_extra,std::ofstream::binary);       //img part  (binary)

                compartment_bout.open(bout_compartment_walkers, std::ofstream::binary); //compartment walkers
            }
            
            std::vector<int> compartment_walkers_intra = std::vector<int>(simulations[0]->dynamicsEngine->list_walkers_intra.size(), 0);
            
            for(unsigned p = 0; p < params.num_proc; p++)
            {
                for (unsigned i = 0 ; i < compartment_walkers_intra.size(); i++){
                    compartment_walkers_intra[i] += simulations[p]->dynamicsEngine->list_walkers_intra[i];
                }
            }
            

            if(params.write_bin){
                std::vector<float> compartment_walkers_intra_float(compartment_walkers_intra.size());
                std::transform(compartment_walkers_intra.begin(), compartment_walkers_intra.end(), compartment_walkers_intra_float.begin(),
                            [](double val) { return static_cast<float>(val); });

                compartment_bout.write(reinterpret_cast<const char *>(compartment_walkers_intra_float.data()), compartment_walkers_intra_float.size() * sizeof(float));
            }
            if(params.write_txt){
                for (unsigned i = 0 ; i < compartment_walkers_intra.size(); i++ ) {
                    compartment_out  << compartment_walkers_intra[i] << " ";
                }
                compartment_out  << std::endl;
            }

            for (unsigned i = 0 ; i < simulations[0]->dataSynth->DWI_intervals_intra.size(); i++)
            {
                std::vector<double> DWI (simulations[0]->dataSynth->DWI_intervals_intra[i].size(), 0);
                std::vector<double> DWIi(simulations[0]->dataSynth->DWIi_intervals_intra[i].size(), 0);

                for (unsigned j = 0 ; j < simulations[0]->dataSynth->DWI_intervals_intra[i].size(); j++ )
                {

                    std::vector<int> phase(3600, 0);
                    for(unsigned p = 0; p < params.num_proc; p++)
                    {
                        DWI[j]+=  simulations[p]->dataSynth->DWI_intervals_intra[i][j];       //real part
                        DWIi[j]+= simulations[p]->dataSynth->DWIi_intervals_intra[i][j];      //img  part

                    }

                    if(params.write_txt){
                        dwi_out_intra  << DWI[j] << " ";
                        dwii_out_intra << DWIi[j] << " ";
                    }
                }

                if(params.write_txt){
                    dwi_out_intra  << std::endl;
                    dwii_out_intra << std::endl;
                }

                if(params.write_bin){
                    std::vector<float> DWI_float(DWI.size());
                    std::transform(DWI.begin(), DWI.end(), DWI_float.begin(),
                                [](double val) { return static_cast<float>(val); });

                    dwi_bout_intra.write(reinterpret_cast<const char *>(DWI_float.data()), DWI_float.size() * sizeof(float));

                    std::vector<float> DWI_floati(DWIi.size());
                    std::transform(DWIi.begin(), DWIi.end(), DWI_floati.begin(),
                                [](double val) { return static_cast<float>(val); });

                    dwii_bout_intra.write(reinterpret_cast<const char *>(DWI_floati.data()), DWI_floati.size() * sizeof(float));

                }

            }
            for (unsigned i = 0 ; i < simulations[0]->dataSynth->DWI_intervals_extra.size(); i++ )
            {
                std::vector<double> DWI (simulations[0]->dataSynth->DWI_intervals_extra[i].size(), 0);
                std::vector<double> DWIi(simulations[0]->dataSynth->DWIi_intervals_extra[i].size(), 0);

                for (unsigned j = 0 ; j < simulations[0]->dataSynth->DWI_intervals_extra[i].size(); j++ )
                {

                    std::vector<int> phase(3600, 0);
                    for(unsigned p = 0; p < params.num_proc; p++)
                    {
                        DWI[j]+=  simulations[p]->dataSynth->DWI_intervals_extra[i][j];       //real part
                        DWIi[j]+= simulations[p]->dataSynth->DWIi_intervals_extra[i][j];      //img  part

                    }

                    if(params.write_txt){
                        dwi_out_extra << DWI[j] << " ";
                        dwii_out_extra << DWIi[j] << " ";
                    }
                }

                if(params.write_txt){
                    dwi_out_extra  << std::endl;
                    dwii_out_extra << std::endl;
                }

                if(params.write_bin){
                    std::vector<float> DWI_float(DWI.size());
                    std::transform(DWI.begin(), DWI.end(), DWI_float.begin(),
                                [](double val) { return static_cast<float>(val); });

                    dwi_bout_extra.write(reinterpret_cast<const char *>(DWI_float.data()), DWI_float.size() * sizeof(float));

                    std::vector<float> DWI_floati(DWIi.size());
                    std::transform(DWIi.begin(), DWIi.end(), DWI_floati.begin(),
                                [](double val) { return static_cast<float>(val); });

                    dwii_bout_extra.write(reinterpret_cast<const char *>(DWI_floati.data()), DWI_floati.size() * sizeof(float));

                }

            }

            if(params.write_txt){
                dwi_out_intra.close();
                dwii_out_intra.close();

                dwi_out_extra.close();
                dwii_out_extra.close();

                compartment_out.close();
            }
            if(params.write_bin){
                dwi_bout_intra.close();
                dwii_bout_intra.close();

                dwi_bout_extra.close();
                dwii_bout_extra.close();

                compartment_bout.close();
            }

            if(params.log_phase_shift){
                phase_out_intra.close();
                phase_out_extra.close();
            }
        }
    
    }

    stuck_count = 0;
    illegal_count = 0;

    for(unsigned i = 0 ; i < simulations.size();i++){
        stuck_count   += simulations[i]->dynamicsEngine->sentinela.stuck_count;
        illegal_count += simulations[i]->dynamicsEngine->sentinela.illegal_count;
        icvf+= simulations[i]->dynamicsEngine->num_simulated_walkers*simulations[i]->dynamicsEngine->icvf;
    }
    icvf/=float(this->total_sim_particles);

    if(params.custom_sampling_area == false and params.voxels_list.size()>0){
        for(auto i = 0; i<3;i++){
           params.min_sampling_area[i]= params.voxels_list[0].first[i];
           params.max_sampling_area[i]= params.voxels_list[0].second[i];
        }
    }

    float sampling_area = 1;

    for (auto i =0; i < 3;i++ ){
        sampling_area*= (params.max_sampling_area[i]-params.min_sampling_area[i]);
    }

    aprox_volumen = (icvf)*sampling_area;



    //====================================================================================================
    //Subdivision results join.
    //====================================================================================================


    if(simulations[0]->dataSynth && this->params.subdivision_flag)
    {
        std::string out_vox_DWI    = params.output_base_name + "_voxels_DWI.txt";
        std::string out_vox_DWIi   = params.output_base_name + "_voxels_DWI_img.txt";

        std::string bout_vox_DWI    = params.output_base_name  + "_voxels_DWI.bfloat";
        std::string bout_vox_DWIi   = params.output_base_name  + "_voxels_DWI_img.bfloat";

        std::ofstream dwi_out, dwii_out, dwi_bout, dwii_bout;

        if(params.write_txt){
            dwi_out.open(out_vox_DWI  ,std::ofstream::out);       //real part
            dwii_out.open(out_vox_DWIi,std::ofstream::out);       //img part
        }

        if(params.write_bin){
            dwi_bout.open(bout_vox_DWI  ,std::ofstream::binary);       //real part (binary)
            dwii_bout.open(bout_vox_DWIi,std::ofstream::binary);       //img part  (binary)
        }


        std::string outPDen   = params.output_base_name  + "_volFractions.txt";
        std::ofstream PDen_out(outPDen,std::ofstream::out);   //Particle Density map

        //## HEADER
        if(params.write_txt){
            // #Num_voxels
            dwi_out  << simulations[0]->dataSynth->subdivisions.size() << endl;
            dwii_out << simulations[0]->dataSynth->subdivisions.size() << endl;
            // #Size_DWI
            dwi_out  << simulations[0]->dataSynth->DWI.size() << endl;
            dwii_out << simulations[0]->dataSynth->DWI.size() << endl;
        }
        //## HEADER (Binary)
        if(params.write_bin){
            // #Num_voxels
            float holder = float(simulations[0]->dataSynth->subdivisions.size());
            dwi_bout.write (reinterpret_cast<char *>(&holder), sizeof(float));
            dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
            // #Size_DWI
            holder = float(simulations[0]->dataSynth->DWI.size());
            dwi_bout.write (reinterpret_cast<char *>(&holder), sizeof(float));
            dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
        }

        for(uint s = 0; s < simulations[0]->dataSynth->subdivisions.size(); s++)
        {
            double density_intra = 0;
            double density_extra = 0;
            double density = 0;
            for(unsigned p = 0; p < params.num_proc; p++){

                density_intra +=  simulations[p]->dataSynth->subdivisions[s].density_intra;
                density_extra +=  simulations[p]->dataSynth->subdivisions[s].density_extra;
                density       +=  simulations[p]->dataSynth->subdivisions[s].density;
            }

            PDen_out << density << endl;
            PDen_out << density_intra << endl;
            PDen_out << density_extra << endl;


            for (unsigned i = 0 ; i < simulations[0]->dataSynth->DWI.size(); i++ )
            {
                double sub_DWI   = 0;
                double sub_DWIi  = 0;

                for(unsigned p = 0; p < params.num_proc; p++)
                {
                    sub_DWI +=  simulations[p]->dataSynth->sub_DWI[s][i];       //real part
                    sub_DWIi+=  simulations[p]->dataSynth->sub_DWIi[s][i];      //img  part
                }


                if(params.write_txt){
                    dwi_out  << sub_DWI  << std::endl;
                    dwii_out << sub_DWIi << std::endl;
                }
                if(params.write_bin){
                    float holder = float(sub_DWI);
                    dwi_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                    holder = float(sub_DWIi);
                    dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                }
            }
        }

        if(params.write_txt){
            dwi_out.close();
            dwii_out.close();
        }
        if(params.write_bin){
            dwi_bout.close();
            dwii_bout.close();
        }

        PDen_out.close();
    }


    //====================================================================================================
    // NEW Join propagator join solution
    //====================================================================================================


    if (this->params.log_propagator){
        uint num_sims = uint(simulations.size());
        for (uint i =0; i< num_sims; i++){

            for (uint j = 0; j < this->simulations[i]->dynamicsEngine->propagator.num_times; j++){

                for (uint k = 0; k < this->simulations[i]->dynamicsEngine->propagator.num_dirs; k++){
                    if(i == 0){
                        this->simulations[0]->dynamicsEngine->propagator.propagator_log[j][k]/=float(num_sims);
                        }
                    else{
                        this->simulations[0]->dynamicsEngine->propagator.propagator_log[j][k]+= this->simulations[i]->dynamicsEngine->propagator.propagator_log[j][k]/float(num_sims);
                    }
                }
            }
        }
    }


    std::string prop_file_name = params.output_base_name + "_propagator.txt";

    if (this->params.log_propagator)
        this->simulations[0]->dynamicsEngine->writePropagator(prop_file_name);
}

void ParallelMCSimulation::specialInitializations()
{
    if(params.gamma_packing == true){

        // Statistic
        GammaDistribution stat_dist(params.packing_num_obstacles, params.gamma_packing_alpha, params.gamma_packing_beta);
        stat_dist.displayDistribution();
        std::vector<double> radiis_list(stat_dist.createRadiiList());
        
        string message = "Initialializing Gamma distribution\n";
        
        for (unsigned i=0;i<params.gamma_packing_alpha.size();++i){
            message+="("+ std::to_string(params.gamma_packing_alpha[i]) + ","
                + std::to_string(params.gamma_packing_beta[i]) + ").\n";
        }
        
        if (params.packing_cyl==true){

            message += " Cylinders \n";
            SimErrno::info(message,cout);        
            
            // Obstacle
            CylinderDistribution packing_dist(params.packing_icvf, params.min_limits, params.max_limits, radiis_list);
            packing_dist.createSubstrate();


            string file = params.output_base_name + "_gamma_packing_cylinder_list.txt";
            params.cylinders_files.push_back(file);
            ofstream out(file);
            packing_dist.printSubstrate(out);
            

            params.max_limits = packing_dist.max_limits;
            params.min_limits = packing_dist.min_limits;

            if(params.voxels_list.size()<=0){
                pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(params.min_limits,params.max_limits);
                params.voxels_list.push_back(voxel_);
            }
            else{
                params.voxels_list[0].first =  params.min_limits;
                params.voxels_list[0].second = params.max_limits;
            }
            SimErrno::info("Done.\n",cout);
        }

        if (params.packing_s==true){
            message += " Spheres \n";
            SimErrno::info(message,cout);      

            // Obstacle
            SphereDistribution packing_dist(params.packing_icvf, params.min_limits, params.max_limits, radiis_list);
            packing_dist.createSubstrate();

            string file = params.output_base_name + "_gamma_packing_sphere_list.txt";
            ofstream out(file);
            packing_dist.printSubstrate(out);            
            params.spheres_files.push_back(file);
            
            params.max_limits = packing_dist.max_limits;
            params.min_limits = packing_dist.min_limits;

            if(params.voxels_list.size()<=0){
                pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(params.min_limits,params.max_limits);
                params.voxels_list.push_back(voxel_);
            }
            else{
                params.voxels_list[0].first =  params.min_limits;
                params.voxels_list[0].second = params.max_limits;
            }

            SimErrno::info("Done.\n",cout);

        } 

    }

    if(params.uniform_packing == true){
        // Statistic
        UniformDistribution stat_dist(params.packing_num_obstacles,params.uniform_packing_radii);
        stat_dist.displayDistribution();
        std::vector<double> radiis_list(stat_dist.createRadiiList());

        string message = "Initialializing Uniform packing of spheres \n";
        for (unsigned i=0;i<params.uniform_packing_radii.size();++i){
            message+="("+ std::to_string(params.uniform_packing_radii[i]) + ").\n";
        }

        if (params.packing_s==true){
            message += " Spheres \n";
            SimErrno::info(message,cout);  

            // Obstacle
            SphereDistribution packing_dist(params.packing_icvf, params.min_limits, params.max_limits, radiis_list);
            packing_dist.createSubstrate();

                        
            string file = params.output_base_name + "_uniform_packing_sphere_list.txt";
            params.spheres_files.push_back(file);
            ofstream out(file);
            packing_dist.printSubstrate(out);

            params.max_limits = packing_dist.max_limits;
            params.min_limits = packing_dist.min_limits;

            if(params.voxels_list.size()<=0){
                pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(params.min_limits,params.max_limits);
                params.voxels_list.push_back(voxel_);
            }
            else{
                params.voxels_list[0].first =  params.min_limits;
                params.voxels_list[0].second = params.max_limits;
            }

            SimErrno::info("Done.\n",cout);
        }

    }
    if(params.gaussian_packing == true){

        // Statistic
        GaussianDistribution stat_dist(params.packing_num_obstacles,params.gaussian_packing_mean, params.gaussian_packing_std);
        stat_dist.displayDistribution();
        std::vector<double> radiis_list(stat_dist.createRadiiList());
        

        string message = "Initialializing Gaussian packing of spheres \n";
        for (unsigned i=0;i<params.gaussian_packing_mean.size();++i){
            message+="("+ std::to_string(params.gaussian_packing_mean[i]) + "," + std::to_string(params.gaussian_packing_std[i])+ ").\n";
        }
        if (params.packing_s==true){

            message += " Spheres \n";
            SimErrno::info(message,cout);  
         
            // Obstacle
            SphereDistribution packing_dist(params.packing_icvf, params.min_limits, params.max_limits, radiis_list);
            packing_dist.createSubstrate();
                         
            string file = params.output_base_name + "_gaussian_packing_sphere_list.txt";
            params.spheres_files.push_back(file);
            ofstream out(file);
            packing_dist.printSubstrate(out);

            params.max_limits = packing_dist.max_limits;
            params.min_limits = packing_dist.min_limits;

            if(params.voxels_list.size()<=0){
                pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(params.min_limits,params.max_limits);
                params.voxels_list.push_back(voxel_);
            }
            else{
                params.voxels_list[0].first =  params.min_limits;
                params.voxels_list[0].second = params.max_limits;
            }

            SimErrno::info("Done.\n",cout);
        }

    }

}








