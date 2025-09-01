#include "trajectory.h"
#include <iomanip>      // std::setprecision
#include "Eigen/src/Core/Matrix.h"
#include "simerrno.h"
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Eigen/src/Core/ArithmeticSequence.h"
#include "Eigen/src/Core/util/IndexedViewHelper.h"
#include "Eigen/src/Core/util/Macros.h"
#include "walker.h"
#include <iostream>
#include <string>
#include <vector>
using namespace std;


Trajectory::Trajectory()
{
    isBigEndian = false;
    in          = NULL;
    in_header   = NULL;
    inhit       = NULL;
    in_headerhit= NULL;
    steps_subset = false;
}


Trajectory:: Trajectory(const char* traj_file, const char* hit_file, const char* fullc_file, bool isBigEndian_, std::string io_flag_)
{
    in              = NULL;
    in_header       = NULL;
    inhit           = NULL;
    in_headerhit    = NULL;
    this->trajfile  = traj_file;
    this->hitfile   = hit_file;
    this->fullcfile = fullc_file;
    
    isBigEndian     = isBigEndian_;
    readTrajectoryHeader();
    readHitHeader();
    io_flag         = io_flag_;
    steps_subset    = false;

}

Trajectory::~Trajectory()
{
    // closeTrajReaderFile();
}

void Trajectory::initTrajectory(Parameters params)
{

    dyn_duration = params.sim_duration;
    T            = params.num_steps;
    N            = params.num_walkers;
    write_traj   = params.write_traj;
    write_txt    = params.write_txt;
    write_bin    = params.write_bin;
    write_hit    = params.write_hit;
    write_full_c = params.write_full_c;
    write_every_nth_step = params.write_every_nth_step;
    write_location = params.write_location;

    trajfile     = params.output_base_name;
    hitfile      = params.output_base_name;
    fullcfile    = params.output_base_name;

    if(params.record_pos_times.size() > 0){
        steps_subset = true;
        pos_times = params.record_pos_times;
        std::sort(pos_times.begin(),pos_times.end());
    }
}


void Trajectory::initTrajWriter()
{
    if(write_traj){
        if(write_bin){
            initTrajWriterBinary();
            writeTrajectoryHeaderBinary();
        }

        if(write_txt){
            initTrajWriterText();
            writeTrajectoryHeaderText();
        }
    }

}


void Trajectory::initHitWriter()
{
    if(write_hit){
        initHitWriterBinary();
        writeHitHeaderBinary();
    }
}

void Trajectory::initFullWriter()
{
    if(write_full_c){
        initFullWriterBinary();
    }
}



void Trajectory::initTrajWriterBinary()
{
    if(bout)
        bout.close();
    if(bheaderout)
        bheaderout.close();

    bout.open((trajfile + ".traj").c_str(), std::ofstream::binary);
    bheaderout.open((trajfile + ".bhdr").c_str(), std::ofstream::binary);

    if(!bout || !bheaderout){
        std::cout << "Cannot open " << trajfile.c_str() << std::endl;
        return;
    }
}

void Trajectory::initHitWriterBinary()
{
    if(bouthit)
        bouthit.close();
    if(bheaderouthit)
        bheaderouthit.close();

    bouthit.open((hitfile + ".hit").c_str(), std::ofstream::binary);
    bheaderouthit.open((hitfile + ".bhdrhit").c_str(), std::ofstream::binary);

    if(!bouthit || !bheaderouthit){
        std::cout << "Cannot open " << hitfile.c_str() << std::endl;
        return;
    }
}


void Trajectory::initFullWriterBinary()
{
    if(boutfull_loc)
        boutfull_loc.close();
    if(boutfull_cross)
        boutfull_cross.close();

    boutfull_loc.open((fullcfile + ".fullloc").c_str(), std::ofstream::binary);
    boutfull_cross.open((fullcfile + ".fullcross").c_str(), std::ofstream::binary);
    
    if(!boutfull_loc || !boutfull_cross){
        std::cout << "Cannot open " << fullcfile.c_str() << std::endl;
        return;
    }
}


void Trajectory::initTrajWriterText()
{
    if(tout)
        tout.close();

    if(theaderout)
        theaderout.close();


    tout.open((trajfile + ".traj.txt").c_str(), std::ios::out);
    theaderout.open((trajfile + ".hdr.txt").c_str(), std::ios::out);

    if(!tout || !theaderout){
        //TODO: Error handling
        std::cout << "Cannot open " << (trajfile + ".traj.txt").c_str()<< std::endl;
        return;
    }

    // write header line 
    if (write_location) {
        tout << "walker_id step_number x y z location" << std::endl;
    }
    else {
        tout << "walker_id step_number x y z" << std::endl;
    }
}



void Trajectory::writePosition(Eigen::Vector3d &pos)
{
    cout << "We use the vector function, caution, subsampling not implemented! Location output as well." ;
    cout << "traj"  << write_traj << endl;
    if(write_traj){
        if(write_bin){
            writePositionBinary(pos);
        }
        if(write_txt){
            writePositionText(pos);
        }
    }

}

void Trajectory::writeTrajectoryHeaderBinary()
{
    // We write everything in float to unify the binary format.
    float duration = float(dyn_duration);
    float N_ = N;

    float T_ = T;
    bheaderout.write(reinterpret_cast<char *>(&duration), sizeof(duration));
    bheaderout.write(reinterpret_cast<char *>(&N_), sizeof(N_));

    if(steps_subset)
        T_ = float(pos_times.size());

    bheaderout.write(reinterpret_cast<char *>(&T_), sizeof(float));
}

void Trajectory::writeHitHeaderBinary()
{
    // We write everything in float to unify the binary format.
    float duration = float(dyn_duration);
    float N_ = N;

    float T_ = T;
    bheaderouthit.write(reinterpret_cast<char *>(&duration), sizeof(duration));
    bheaderouthit.write(reinterpret_cast<char *>(&N_), sizeof(N_));

    if(steps_subset)
        T_ = float(pos_times.size());

    bheaderouthit.write(reinterpret_cast<char *>(&T_), sizeof(float));
}

void Trajectory::writeTrajectoryHeaderText()
{
    theaderout << dyn_duration << std::endl;
    theaderout << N << std::endl;

    if(steps_subset)
        theaderout << this->pos_times.size() << std::endl;
    else
        theaderout << T << std::endl;
}

void Trajectory::reWriteHeaderFile(unsigned num_walkers)
{
    if(bheaderout)
        bheaderout.close();

    if(theaderout)
        theaderout.close();

    if(bheaderouthit)
        bheaderouthit.close();

    if(write_traj)
    {
        if(write_txt){
            theaderout.open((trajfile + ".hdr.txt").c_str(), std::ios::out);
            if( !theaderout){
                //TODO: Error handling
                std::cout << "Cannot open header:  " << trajfile.c_str() << std::endl;
                return;
            }

            theaderout << dyn_duration << std::endl;
            theaderout << num_walkers << std::endl;

            if(steps_subset)
                theaderout << this->pos_times.size() << std::endl;
            else
                theaderout << T << std::endl;
        }

        if(write_bin){
            bheaderout.open((trajfile + ".bhdr").c_str(), std::ofstream::binary);
            if( !bheaderout){
                //TODO: Error handling
                std::cout << "Cannot open header: " << trajfile.c_str() << std::endl;
                return;
            }

            // We write everything in float to unify the binary format.
            float duration = float(dyn_duration);
            float N_ = num_walkers;

            float T_ = T;
            bheaderout.write(reinterpret_cast<char *>(&duration), sizeof(duration));
            bheaderout.write(reinterpret_cast<char *>(&N_), sizeof(N_));

            if(steps_subset)
                T_ = float(pos_times.size());

            bheaderout.write(reinterpret_cast<char *>(&T_), sizeof(float));
        }
    }

        if(write_hit)
    {
        bheaderouthit.open((hitfile + ".bhdrhit").c_str(), std::ofstream::binary);
        if( !bheaderouthit){
            //TODO: Error handling
            std::cout << "Cannot open header: " << hitfile.c_str() << std::endl;
            return;
        }

        // We write everything in float to unify the binary format.
        float duration = float(dyn_duration);
        float N_ = num_walkers;

        float T_ = T;
        bheaderouthit.write(reinterpret_cast<char *>(&duration), sizeof(duration));
        bheaderouthit.write(reinterpret_cast<char *>(&N_), sizeof(N_));

        if(steps_subset)
            T_ = float(pos_times.size());

        bheaderouthit.write(reinterpret_cast<char *>(&T_), sizeof(float));
        
    }
}

void Trajectory::writePositionText(Eigen::Vector3d &pos)
{
    tout << std::setprecision(6) << pos[0] << std::endl << pos[1] << std::endl << pos[2] << std::endl;
}

void Trajectory::writePositionBinary(Eigen::Vector3d &pos)

{
    float pos0 = float(pos(0)),pos1 = float(pos(1)),pos2 = float(pos(2));
    bout.write(reinterpret_cast<char *>(&pos0), sizeof(float));
    bout.write(reinterpret_cast<char *>(&pos1), sizeof(float));
    bout.write(reinterpret_cast<char *>(&pos2), sizeof(float));
}


void Trajectory::writePosition(Eigen::Matrix3Xd &pos, Eigen::VectorXi &col_in, Eigen::VectorXi &col_ext, Eigen::VectorXi &cross_in, Eigen::VectorXi &cross_ext)
{

    if(write_traj)
    {
        unsigned number_of_subsampled_columns = int(pos.cols() / write_every_nth_step)+1;
        Eigen::Matrix3Xd pos_subsampled(3, number_of_subsampled_columns);

        for (int c = 0; c < number_of_subsampled_columns; ++c) {
            pos_subsampled.col(c) = pos.col(c * write_every_nth_step);
        }

        // add last position
        if (pos.cols() % write_every_nth_step != 0) {
            pos_subsampled.conservativeResize(Eigen::NoChange, pos_subsampled.cols() + 1);
            pos_subsampled.col(pos_subsampled.cols() - 1) = pos.col(pos.cols() - 1);
        }

        if(write_bin)
            writePositionBinary(pos_subsampled);

        if(write_txt)
            writePositionText(pos_subsampled);
    }

    if(write_hit)
    {
        writePositionHit(col_in, col_ext, cross_in, cross_ext);
    }    
}  

void Trajectory::writePosition(Walker &walker, unsigned &walker_index){
    // Malte's new writePosition for text output that supports location output (intra/extra)
    Eigen::Matrix3Xd& pos = walker.pos_r_log;

    if (write_traj){
        if(write_bin){
            cout << "Only text output implemented for location output." << std::endl;

        }
        else if (write_txt) {
            if (write_every_nth_step>1) {
                unsigned number_of_subsampled_columns = int(pos.cols() / write_every_nth_step)+1;
                Eigen::Matrix3Xd pos_subsampled(3, number_of_subsampled_columns);
                std::vector<int> step_number;

                for (int c = 0; c < number_of_subsampled_columns; ++c) {
                    pos_subsampled.col(c) = pos.col(c * write_every_nth_step);
                    step_number.push_back(write_every_nth_step*c);
                }

                // add last position
                if (pos.cols() % write_every_nth_step != 0) {
                    pos_subsampled.conservativeResize(Eigen::NoChange, pos_subsampled.cols() + 1);
                    pos_subsampled.col(pos_subsampled.cols() - 1) = pos.col(pos.cols() - 1);
                    step_number.push_back(T+1);
                }
                            
                for (int i = 0; i < pos_subsampled.cols(); i++){
                    if (write_location) {
                        std::string this_location;
                        if (walker.location == Walker::intra){
                            this_location="intra";
                        }
                        else if (walker.location == Walker::extra){
                            this_location="extra";
                        }
                        else if (walker.location == Walker::unknown){
                            this_location="unknown";
                        }
                        else {
                            this_location="not_implemented";
                        }

                        tout << std::setprecision(6) 
                            << walker_index << " " 
                            << step_number[i] << " "
                            << pos(0,i) << " " << pos(1,i) << " " 
                            << pos(2,i) << " " 
                            << this_location << std::endl;
                    }
                    else {
                        tout << std::setprecision(6) 
                             << walker_index << " " 
                             << step_number[i] << " "
                             << pos(0,i) << " " 
                             << pos(1,i) << " " 
                             << pos(2,i) << std::endl;
                }
            }
            }
            




        }
    }
}
    

void Trajectory::writeFullCollision(Eigen::Vector3d &col_point, int &cross, int &loc, unsigned &t, unsigned &id_)
{
    if(write_full_c)
    {
        float pos0 = float(col_point(0)),pos1 = float(col_point(1)),pos2 = float(col_point(2)) , t0 = float(t), id_0 = float(id_);
                boutfull_loc.write(reinterpret_cast<char *>(&pos0), sizeof(float));
                boutfull_loc.write(reinterpret_cast<char *>(&pos1), sizeof(float));
                boutfull_loc.write(reinterpret_cast<char *>(&pos2), sizeof(float));
                boutfull_loc.write(reinterpret_cast<char *>(&t0), sizeof(float));
                boutfull_loc.write(reinterpret_cast<char *>(&id_0), sizeof(float));

        int cross0 = static_cast<int8_t>(cross), loc0 = static_cast<int8_t>(loc) ;
            boutfull_cross.write(reinterpret_cast<char *>(&cross0), sizeof(int8_t));
            boutfull_cross.write(reinterpret_cast<char *>(&loc0), sizeof(int8_t));

    }
}


void Trajectory::writePositionText(Eigen::Matrix3Xd &pos)
{
    cout << "Warning, buggy" << std::endl;
    if(steps_subset)
    {
        unsigned index = 0;
        for(unsigned i = 0; i < pos.cols(); i++ ){
            if(i == pos_times[index]){
                tout << std::setprecision(6) << pos(0,i) << std::endl << pos(1,i) << std::endl << pos(2,i) << std::endl << std::endl;
                index++;  
            }
            else
            {
                for(unsigned i = 0; i < pos.cols(); i++ ){
                    tout << std::setprecision(6) << pos(0,i) << std::endl << pos(1,i) << std::endl << pos(2,i) << std::endl << std::endl;
                }
            }
        }
    }
    else
    {
        for(unsigned  i = 0; i < pos.cols(); i++ ){      
            float pos0 = float(pos(0,i)),pos1 = float(pos(1,i)),pos2 = float(pos(2,i));
            bout.write(reinterpret_cast<char *>(&pos0), sizeof(float));
            bout.write(reinterpret_cast<char *>(&pos1), sizeof(float));
        }
    }

}

// Prefer a const reference and Eigen::Index for indexing
void Trajectory::writePositionBinary(const Eigen::Matrix3Xd& pos)
{
    assert(pos.rows() == 3);

    if (pos.cols() == 0) return;  // nothing to write

    auto write_col = [&](Eigen::Index c)
    {
        // c must be valid
        assert(c >= 0 && c < pos.cols());
        const Eigen::Vector3f p = pos.col(c).cast<float>();
        bout.write(reinterpret_cast<const char*>(p.data()), 3 * sizeof(float));
    };

    if (!steps_subset) {
        // ✅ NEVER use T here — write exactly what's in the matrix
        for (Eigen::Index c = 0; c < pos.cols(); ++c)
            write_col(c);
        return;
    }

    // steps_subset == true
    if (pos_times.empty()) return;

    // Optional sanity: ensure non-negative and (usually) sorted
    // assert(std::is_sorted(pos_times.begin(), pos_times.end()));

    // Decide whether 'pos' stores all timesteps (full) or only the kept ones (compressed)
    const bool full_layout = (pos.cols() == static_cast<Eigen::Index>(T + 1));

    if (full_layout) {
        // Use actual timestep as column
        for (size_t k = 0; k < pos_times.size(); ++k) {
            const long t = static_cast<long>(pos_times[k]);  // 0-based expected
            if (t < 0 || t >= pos.cols()) {
                // out-of-range timestep in config; skip or throw
                // std::cerr << "Warning: timestep " << t << " >= pos.cols()=" << pos.cols() << '\n';
                continue;
            }
            write_col(static_cast<Eigen::Index>(t));
        }
    } else {
        // Compressed layout: kth kept time is at column k
        const Eigen::Index want = static_cast<Eigen::Index>(pos_times.size());
        if (pos.cols() != want) {
            // Layout mismatch; write the common part without crashing
            const Eigen::Index n = std::min(pos.cols(), want);
            // std::cerr << "Warning: pos.cols()=" << pos.cols() << " != pos_times.size()=" << want << '\n';
            for (Eigen::Index c = 0; c < n; ++c) write_col(c);
        } else {
            for (Eigen::Index c = 0; c < pos.cols(); ++c) write_col(c);
        }
    }
}



void Trajectory::writePositionHit(Eigen::VectorXi &col_in, Eigen::VectorXi &col_ext, Eigen::VectorXi &cross_in, Eigen::VectorXi &cross_ext)
{
    if(steps_subset)
    {
        unsigned index = 0;
        for(unsigned i = 0; i < T+1; i++ )
            if(i == pos_times[index]){

                int col0 = static_cast<int8_t>(col_in(i)), col1 = static_cast<int8_t>(col_ext(i)), cross0 = static_cast<int8_t>(cross_in(i)), cross1 = static_cast<int8_t>(cross_ext(i));
                bouthit.write(reinterpret_cast<char *>(&col0), sizeof(int8_t));
                bouthit.write(reinterpret_cast<char *>(&col1), sizeof(int8_t));
                bouthit.write(reinterpret_cast<char *>(&cross0), sizeof(int8_t));
                bouthit.write(reinterpret_cast<char *>(&cross1), sizeof(int8_t));
                index++;                        // Update the index

                if(index >= pos_times.size()){
                    break;
                }
            }
    }
    else
    {   
        /*
        for(unsigned  i = 0; i < T+1; i++ ){      
            int col0 = static_cast<int8_t>(col_in(i)), col1 = static_cast<int8_t>(col_ext(i)), cross0 = static_cast<int8_t>(cross_in(i)), cross1 = static_cast<int8_t>(cross_ext(i));
                bouthit.write(reinterpret_cast<char *>(&col0), sizeof(int8_t));
                bouthit.write(reinterpret_cast<char *>(&col1), sizeof(int8_t));
                bouthit.write(reinterpret_cast<char *>(&cross0), sizeof(int8_t));
                bouthit.write(reinterpret_cast<char *>(&cross1), sizeof(int8_t));
            }

        */
       // Short version - In combination with full colision
        for(unsigned  i = 0; i < T+1; i++ ){      
            int col0 = static_cast<int8_t>(col_in(i) +col_ext(i)), cross0 = static_cast<int8_t>(cross_in(i)+cross_ext(i));
                bouthit.write(reinterpret_cast<char *>(&col0), sizeof(int8_t));
                bouthit.write(reinterpret_cast<char *>(&cross0), sizeof(int8_t));
            }
                
    }
}


void Trajectory::initTrajReaderFile()
{
    openTrajReaderFile();
    //fseek ( in , 24 , SEEK_SET );
}

void Trajectory::initHitReaderFile()
{
    openHitReaderFile();
    //fseek ( in , 24 , SEEK_SET );
}

void Trajectory::openTrajReaderFile()
{

    if(io_flag == ""){
        io_flag = "rb";
    }
    //closeTrajReaderFile();

    in        = fopen(trajfile.c_str()  ,io_flag.c_str());
    in_header = fopen(headerfile.c_str(),io_flag.c_str());

    if(!in){
        SimErrno::error("ERROR opening:  " + trajfile + " flag:" + io_flag,std::cout);
        assert(0);
        return;
    }

    if( !in_header){
        SimErrno::error("ERROR opening:  " + headerfile + " flag: " + io_flag,std::cout);
        assert(0);
        return;
    }
}


void Trajectory::openHitReaderFile()
{

    if(io_flag == ""){
        io_flag = "rb";
    }
    //closeTrajReaderFile();

    inhit        = fopen(hitfile.c_str()  ,io_flag.c_str());
    in_headerhit = fopen(hitfile.c_str(),io_flag.c_str());

    if(!inhit){
        SimErrno::error("ERROR opening:  " + hitfile + " flag:" + io_flag,std::cout);
        assert(0);
        return;
    }

    if( !in_headerhit){
        SimErrno::error("ERROR opening:  " + headerfilehit + " flag: " + io_flag,std::cout);
        assert(0);
        return;
    }
}


void Trajectory::closeTrajReaderFile()
{
    if(in != NULL)
        fclose(in);

    if(in_header != NULL)
        fclose(in_header);
}

void Trajectory::closeHitReaderFile()
{
    if(inhit != NULL)
        fclose(inhit);

    if(in_headerhit != NULL)
        fclose(in_headerhit);
}

void Trajectory::setTrajFile(std::string trajfile_)
{

    trajfile = trajfile_  + ".traj";
    headerfile = trajfile_+ ".bhdr";
    readTrajectoryHeader();
}


void Trajectory::setHitFile(std::string hitfile_)
{
    hitfile = hitfile_  + ".hit";
    headerfilehit = hitfile_+ ".bhdrhit";
    readHitHeader();
}

void Trajectory::readTrajectoryHeader()
{

    ifstream myfile (this->headerfile, ios::binary);

    if(!myfile){
        SimErrno::error("Input trajfile header not found!",std::cout);
        assert(0);
    }

    streampos begin,end;
    begin = myfile.tellg();
    myfile.seekg (0, ios::end);
    end 	= myfile.tellg();

    unsigned size_ = end - begin;

    if(size_!= 12){
        SimErrno::error("Corrupted header file",std::cout);
        //cout << "[WARNING] Corrupted header file" << endl;
        myfile.close();
        assert(0);
        return;
    }

    char * memblock = new char[12];

    myfile.close();
    myfile.open(this->headerfile, ios::binary);

    myfile.read(memblock,3*sizeof(float));


    float* float_values = (float*)memblock; //reinterpret as float

    this->N 			= uint(float_values[1]);
    this->T 			= uint(float_values[2]);
    this->dyn_duration 	= double(float_values[0]);

    myfile.close();
    delete[] memblock;

}


void Trajectory::readHitHeader()
{

    ifstream myfile (this->headerfilehit, ios::binary);

    if(!myfile){
        SimErrno::error("Input trajfile header not found!",std::cout);
        assert(0);
    }

    streampos begin,end;
    begin = myfile.tellg();
    myfile.seekg (0, ios::end);
    end 	= myfile.tellg();

    unsigned size_ = end - begin;

    if(size_!= 12){
        SimErrno::error("Corrupted header file",std::cout);
        //cout << "[WARNING] Corrupted header file" << endl;
        myfile.close();
        assert(0);
        return;
    }

    char * memblock = new char[12];

    myfile.close();
    myfile.open(this->headerfilehit, ios::binary);

    myfile.read(memblock,3*sizeof(float));


    float* float_values = (float*)memblock; //reinterpret as float

    this->N 			= uint(float_values[1]);
    this->T 			= uint(float_values[2]);
    this->dyn_duration 	= double(float_values[0]);

    myfile.close();
    delete[] memblock;

}


void Trajectory::readCurrentWalkersTrajectory(Eigen::Matrix3Xd &steps_log)
{
    float x,y,z;

    for (unsigned  t = 0; t < T+1 ; t++){

        assert(fread(&x,4,1,in));
        assert(fread(&y,4,1,in));
        assert(fread(&z,4,1,in));

        if(isBigEndian){
            swapBE2SE2(&x, sizeof(x));
            swapBE2SE2(&y, sizeof(y));
            swapBE2SE2(&z, sizeof(z));
        }

        steps_log(0,t) = x;
        steps_log(1,t) = y;
        steps_log(2,t) = z;
    }
}



//swap between big endian to little endian.
void Trajectory::swapBE2SE2(void *source, int size)
{
    typedef unsigned char TwoBytes[2];
    typedef unsigned char FourBytes[4];
    typedef unsigned char EightBytes[8];

    unsigned char temp;

    if(size == 2)
    {
        TwoBytes *src = (TwoBytes*) source;
        temp = (*src)[0];
        (*src)[0] = (*src)[1];
        (*src)[1] = temp;

        return;
    }

    if(size == 4)
    {
        FourBytes *src = (FourBytes *)source;
        temp = (*src)[0];
        (*src)[0] = (*src)[3];
        (*src)[3] = temp;

        temp = (*src)[1];
        (*src)[1] = (*src)[2];
        (*src)[2] = temp;

        return;
    }

    if(size == 8)
    {
        EightBytes *src = (EightBytes *)source;
        temp = (*src)[0];
        (*src)[0] = (*src)[7];
        (*src)[7] = temp;

        temp = (*src)[1];
        (*src)[1] = (*src)[6];
        (*src)[6] = temp;

        temp = (*src)[2];
        (*src)[2] = (*src)[5];
        (*src)[5] = temp;

        temp = (*src)[3];
        (*src)[3] = (*src)[4];
        (*src)[4] = temp;

        return;
    }

}