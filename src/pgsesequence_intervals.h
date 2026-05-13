//!  PGSE Sequence Class  =============================================================/
/*!
  Derived Class.
  Implementation of the PGSE protocol

  \date   May 2016
  \author Jonathan Rafael
  \version 0.1.0
*=====================================================================================*/

#ifndef PGSEINTSEQUENCE_H
#define PGSEINTSEQUENCE_H

#include <vector>
#include <iostream>
#include "trajectory.h"
#include <Eigen/Core>
#include "simulablesequence.h"
#include "scheme.h"
typedef unsigned long ulong;


/*! \class  ParallelMCSimulation
 *  \brief  Implementation of the PGSE protocol
 */
class PGSESequence_Intervals: public SimulableSequence{
public:

    double TE;           /*!< Time Echo.                                        */

    int T;               /*!< num bins (time steps)                             */

    std::vector<std::vector<double>> t_intervals;    

    std::vector<int> nbr_intervals; /*!< number of intervals in the scheme file             */

    double dyn_duration; /*!< simulation duration (miliseconds)                 */

    std::vector< std::vector<double> > scheme;  /*!< Scheme file values         */

    Trajectory trajectory;    /*!< If the signal is computed from a .trajfile   */

    std::vector<std::vector<std::vector<double>>> sub_DWI_intervals_intra; /*!< auxiliar DWI for signal computations.  */

    std::vector<std::vector<std::vector<double>>> sub_DWI_intervals_extra; /*!< auxiliar DWI for signal computations.  */

    std::vector<std::vector<std::vector<double>>> sub_DWIi_intervals_intra; /*!< auxiliar DWI for signal computations.  */

    std::vector<std::vector<std::vector<double>>> sub_DWIi_intervals_extra; /*!< auxiliar DWI for signal computations.  */

    bool update_DWI;

    int interval_nbr;


    //constructors

    /**
     * @brief Default constructor, set default NULL values. Not to be used.
     */
    PGSESequence_Intervals();
    /**
     * @brief Main constructor. Takes a pre-loaded Scheme file.
     */
    PGSESequence_Intervals(Scheme scheme_);
    /**
     * @brief Main constructor. Takes a pre-loaded Scheme file and a traj file name.
     *        if this argument is passed a traj file is should be written.
     */
    PGSESequence_Intervals(Scheme scheme_,const char* traj_file_name);
    /**
     * @brief Constructor. Takes a the scheme file name to be loaded.
     */
    PGSESequence_Intervals(const char* scheme_file_name);
    /**
     * @brief Constructor. Takes a scheme file name to be loaded and atraj file name.
     *        if this argument is passed a traj file is should be written.
     */
    PGSESequence_Intervals(const char* scheme_file_name,const char* traj_file_name);

    /**
     * @brief Destuctor. Does nothing.
     */
    ~PGSESequence_Intervals();

    /**
     * @brief For using w/o the adt array
     */
    void getGradImpulse(int i, double t, double tLast, Eigen::Vector3d &Gdt) override;


    /**
     * @brief Analytical defined b-value
     */
    double getbValue(unsigned);

    /**
     * @brief Expected free Decay
     */
    double getFreeDecay(unsigned i,double D);


    //TODO: to implement
    /**
     * @brief \warning not implemented yet.
     */
    double getNumericalbValue(unsigned);

    /**
     * @brief Computes de DW signal from a trajfile
     */
    void getDWISignal();

    /**
     * @brief reads the scheme files
     */
    void readSchemeFile();

    /**
     * @param i: updated walker
     */
    virtual void update_phase_shift(double dt,double dt_last, const Walker& walker) override;

    /**
     * @brief Updates the phase shift using the full stored trajectory
     */
    virtual void update_phase_shift(double time_step, const Eigen::Matrix3Xd& trajectory) override;

    /**
     * @brief Updates the DWI signal using the cumulated phase shift
     */
    virtual void update_DWI_signal(const Walker& walker) override;

    /**
     * @brief computes de signal value and sign in a certain time step.
     */
    double get_adt(int grad_index, double t, double tLast);

    /**
     * @brief prints the array adt in the format ().
     */
    double print_adt_and_dt(int grad_index, double t, double tLast);

    virtual void setNumberOfSteps(unsigned T);

    virtual void computeDynamicTimeSteps();

    void SetTimingsIntervals();


private:
    virtual void readSchemeParameters(Scheme scheme_);



};
#endif // PGSEINTSEQUENCE_H
