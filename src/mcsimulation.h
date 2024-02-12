
//███╗   ███╗ ██████╗    ██╗██████╗  ██████╗    ███████╗██╗███╗   ███╗██╗   ██╗██╗      █████╗ ████████╗ ██████╗ ██████╗
//████╗ ████║██╔════╝   ██╔╝██╔══██╗██╔════╝    ██╔════╝██║████╗ ████║██║   ██║██║     ██╔══██╗╚══██╔══╝██╔═══██╗██╔══██╗
//██╔████╔██║██║       ██╔╝ ██║  ██║██║         ███████╗██║██╔████╔██║██║   ██║██║     ███████║   ██║   ██║   ██║██████╔╝
//██║╚██╔╝██║██║      ██╔╝  ██║  ██║██║         ╚════██║██║██║╚██╔╝██║██║   ██║██║     ██╔══██║   ██║   ██║   ██║██╔══██╗
//██║ ╚═╝ ██║╚██████╗██╔╝   ██████╔╝╚██████╗    ███████║██║██║ ╚═╝ ██║╚██████╔╝███████╗██║  ██║   ██║   ╚██████╔╝██║  ██║
//╚═╝     ╚═╝ ╚═════╝╚═╝    ╚═════╝  ╚═════╝    ╚══════╝╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝


//!  Aplication Main Class ======================================================================================/
/*!
*   \details   Main implementation class. Incorporates the particle's dynamics and the data synthesis.
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42.14_wf
*===============================================================================================================*/

#ifndef MCSIMULATION_H
#define MCSIMULATION_H

#include "dynamicsSimulation.h"
#include "scheme.h"


/*! \class MCSimulation
 * \brief  Main implementation class. Incorporates the particle's dynamics and the data synthesis.
 */
class MCSimulation
{
public:
    static int count;                   /*!< count of                                                       */

    int id;                             /*!< Unique id of the simulation                                    */

    DynamicsSimulation* dynamicsEngine; /*!< Instance for the particle dynamics                             */

    Scheme scheme;                      /*!< Scheme file, only PGSE in camino format is supported in 0.2    */

    Parameters params;                  /*!< Parameters instance1 \see :Parameters:                         */

    SimulableSequence* dataSynth;       /*!< Simuleable sequence instance, only PGSE is supported in 0.2    */


    /*! \fn  MCSimulation.
     *  \brief  Default constructor. Intialize everything with 0's and NULL states, object indexes are set to -1.
     */
    MCSimulation();


    /*! \fn  MCSimulation
     *  \param config_file .conf file name with the full set of experiments parameters.
     *  \brief Main constructor.
     */
    MCSimulation(std::string config_file);

    /*! \fn  MCSimulation
     *  \param params_ preloaded simulation parameters.
     *  \brief Secondary constructor.
     */
    MCSimulation(Parameters &params_);

    /*! \fn  ~MCSimulation.
     *  \brief Main destructor. Frees dynamicly allocated memory instances.
     */
    ~MCSimulation();

    /*! \fn  startSimulation
     *  \brief Warp function. Calls the dynamicEngine's native DynamicsSimulation::startSimulation function. \see :DynamicsSimulation:.
     */
    void startSimulation();

    /*! \fn getExpectedFreeeDecay(
     *  \param i, index of the shell/direction on the scheme file.
     *  \brief return the expected free decay based on the parameters of the shell \a i in the scheme file.
     */
    double getExpectedFreeeDecay(unsigned i);

    /*!
     *  Adds all the obstacles defined in the confiuration files.
     */
    void iniObstacles();

private:

    void addCylindersObstaclesFromFiles();

    void addAxonsObstaclesFromFiles();

    void addPLYObstaclesFromFiles();

    void addVoxels();

    void addCylindersConfigurations();

    void addExtraObstacles();

    void addSpheresObstaclesFromFiles();

};

#endif // MCSIMULATION_H
