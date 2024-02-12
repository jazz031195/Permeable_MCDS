//!  Auxiliary class. Handles i/o operation of walker trayectories. ============================/
/*!
*   \author Jonathan Rafael
*   \date   July 2016
*   \version   0.2
*=============================================================================================*/
#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <Eigen/Sparse>
#include <parameters.h>


class Trajectory {
public:
    std::string trajfile;               /*!< trajfile name                                          */
    std::string headerfile;             /*!< header name                                            */
    std::string hitfile;                /*!< hitfile name                                           */
    std::string fullcfile;              /*!< fullfile name                                           */
    
    std::string headerfilehit;          /*!< header name                                            */
    FILE* in, *in_header;               /*!< Files to be written using the previous names           */
    FILE* inhit, *in_headerhit;         /*!< Files to be written using the previous names           */

    /*!< binary out, text out, binary header, text header */
    std::ofstream bout,tout,bheaderout,theaderout, bouthit, bheaderouthit, boutfull_loc, boutfull_cross;


    unsigned N,T;                         /*!< number of walkers, total time;                       */
    //dynamic duration.
    double dyn_duration;
    std::string io_flag;

    std::vector<unsigned> pos_times;    /*!< Times indexes when to save the particle positions.     */

    bool isBigEndian;                   /*!< flag if the format is big endian                       */
    bool write_traj;                    /*!< flag if we want to write a traj file                   */
    bool write_hit;                     /*!< flag if we want to write a hit file                    */
    bool write_txt;                     /*!< flag if we want to write a text traj file              */
    bool write_bin;                     /*!< flag if we want to write a binary traj file            */

    bool write_full_c;

    bool steps_subset;                  /*!< true if the steps are no uniform                       */

    /*!
     * \brief   Main constructor, Initialice everythin to default.
     */
    Trajectory();

    /*!
     * \brief   Contructor , Initialice everythin by parameters.
     */
    Trajectory(const char* traj_file, const char* hit_file, const char* fullc_file, bool isBigEndian_ = true, std::string io_flag_= "rb");

    /*!
     * \brief   Destructor, close files and fstreams.
     */
    ~Trajectory();

    /*!
     * \brief   Initialice the output files if any.
     */
    void initTrajectory(Parameters params);

    /*!
     * \brief   Setd the traj file operations
     */
    void setTrajFile(std::string);

    /*!
     * \brief   Setd the traj file operations
     */
    void setHitFile(std::string);


    //@{
    /** Read operations */
    void closeTrajReaderFile();
    void openTrajReaderFile();
    void initTrajReaderFile();
    void readTrajectoryHeader() ;
    void readCurrentWalkersTrajectory(Eigen::Matrix3Xd&);


    void closeHitReaderFile();
    void openHitReaderFile();
    void initHitReaderFile();
    void readHitHeader() ;
    //@}

    //@{
    /** Write operations */
    void initTrajWriter();
    void initTrajWriterBinary();
    void initTrajWriterText();

    void writeTrajectoryHeaderBinary();
    void writeTrajectoryHeaderText();


    void initHitWriter();
    void initHitWriterBinary();
    
    void writeHitHeaderBinary();

    void initFullWriter();
    void initFullWriterBinary();


    void reWriteHeaderFile(unsigned num_walkers);

    void writePosition(Eigen::Vector3d&);
    void writePositionText(Eigen::Vector3d&);
    void writePositionBinary(Eigen::Vector3d&);

    void writePosition(Eigen::Matrix3Xd&, Eigen::VectorXi&, Eigen::VectorXi&, Eigen::VectorXi&, Eigen::VectorXi&);
    void writePositionText(Eigen::Matrix3Xd&);
    void writePositionBinary(Eigen::Matrix3Xd&);
    void writePositionHit(Eigen::VectorXi&, Eigen::VectorXi&, Eigen::VectorXi&, Eigen::VectorXi&);

    void writeFullCollision(Eigen::Vector3d&, int&, int&, unsigned&, unsigned&);

    //@}

private:
    void swapBE2SE2(void *source, int size);
};

#endif // TRAJECTORY_H
