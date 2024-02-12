//!  Spin Final class =============================================================/
/*!
  Basic unit of the diffusion process.

  \author Jonathan Rafael
  \date November 2016
  \version 0.2
====================================================================================*/

#ifndef WALKER_H
#define WALKER_H

#include "Eigen/Core"
#include <vector>
#include <deque>
#include "collisionsphere.h"
#include <iostream>

/*! \class Walker
 *  \brief Alias for a particle. Basic unit on the simulation process. Saves all the necessary information to perform
 *         the particles dynamics.
 */
class Walker {
public:
    //! An enum.
    /*! All the possibles states that a walker can be in a given step.
        The next step is perform according to this state
    */
    enum state {on_object, on_edge, on_vertex, on_voxel, free, bouncing};

    //! An enum.
    /*! Possible location of the walker inside the voxel.
        Checks illegal crossings of the barrier (border, lol)
    */
    enum RelativeLocation{unknown,intra,extra};

    Eigen::Vector3d pos_r;                                          /*!< Real walker position for collision, r stands for real                  */

    Eigen::Vector3d pos_v;                                          /*!< Walker current position                                                */

    Eigen::Vector3d last_pos_r;                                     /*!< Walker voxel last position                                             */

    Eigen::Vector3d last_pos_v;                                     /*!< Walker real last position                                              */

    Eigen::Vector3d ini_pos;                                        /*!< Walker intital position                                                */

    Eigen::Vector3d next_direction;                                 /*!< Auxiliar vector for special states cases, decides the next direction   */

    Eigen::Matrix3Xd pos_r_log;                                     /*!< log of the real spin position, used to compute the phase shift         */

    Eigen::Matrix3Xd pos_v_log;                                     /*!< log of the voxel position, used for collision location and bouncing    */

    int in_obj_index;                                               /*!< Auxiliar index to save if the walker was inside a convex object        */

    CylinderCollisionSphere collision_sphere_cylinders;             /*!< Collision sphere for collition against cylidners                       */

    AxonCollisionSphere collision_sphere_axons;                     /*!< Collision sphere for collition against axons                      */

    PLYCollisionSphere collision_sphere_ply;                        /*!< Collision sphere for collition against PLY meshes                      */

    SphereCollisionSphere collision_sphere_spheres;                 /*!< Collision sphere for collition against spheres                         */

    Eigen::Vector3d initial_sphere_pos_v;                           /*!< Saves the intial positioon of the walker inside the collition sphere   */

    unsigned steps_count;                                           /*!< Counts the number of steps (including bouncings) made.                 */

    state status;                                                   /*!< state memeber                                                          */

    RelativeLocation initial_location, location, previous_location;                    /*!< location on the substrate (if known)                                   */

    Eigen::VectorXi colision_in_log, colision_ext_log;                                   /*!< Vector of colision for logging                                         */

    Eigen::VectorXi crossing_in_log, crossing_ext_log;                                   /*!< Vector of crossing for logging                                         */

    unsigned colision_in, colision_ext;                                              /*!< Retains the number of hit per step                                      */ 

    unsigned crossing_in, crossing_ext;                                              /*!< Retains the number of crossing per step                                */

    int intra_extra_consensus;                                      /*!< intra o extra position by face collision consensus w/r the normal*/

    unsigned intra_coll_count;                                      /*!< counter of collision in the ïntra-side w/r the normal*/

    unsigned extra_coll_count;                                      /*!< counter of collision in the extra-side w/r the normal*/

    unsigned int index;                                             /*!< Walker identifier (id)*/

    unsigned int rejection_count;                                   /*!< counter of the rejected directions in a single time-step*/

    float steps_per_second;                                         /*!< Particles steps per second speeed.*/

    int in_ax_index;                                                /*!< Index of axon in which walker is inside (-1 if is in extracellular space)*/

    int in_sph_index;                                               /*!< Index of sphere in which walker is inside (-1 if is in extracellular space)*/
    
    bool is_allowed_to_cross;                                       /*!< Is allowed to cross membrane, not illegal crossing */
    
    std::vector<Eigen::Vector3d> normals;                                         /*!< Normal vector for the collisions against voxel boundaries */
    //! Default constructor.
    /*! Set all variables to cero.*/
    Walker();

    //! Default destructor.
    //!
    /*! Does nothing
    */
    ~Walker() {}

    //! Constructor.
    /*! Initialize the walker position in a random position inside the boundaries
        defined by the limits.
        \param xmin lower x threshold
        \param xmax upper x threshold
        \param ymin lower y threshold
        \param ymax upper y threshold
    */
    Walker(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

    // Get methods
    void  getRealPosition(double &, double &, double &) const;
    void  getRealPosition(Eigen::Vector3d&) const;
    void  getVoxelPosition(double &, double &, double &) const;
    void  getVoxelPosition(Eigen::Vector3d&) const;
    void  getInitialPosition(double &, double &, double &) const;
    void  getInitialPosition(Eigen::Vector3d&) const;
    void  getNextDirection(Eigen::Vector3d&) const;
    unsigned int getIndex() const;

    // Set methods
    void  setRealPosition(const double &, const double &,const double &);
    void  setRealPosition(const Eigen::Vector3d&);
    void  setVoxelPosition(const double &, const double &,const double &);
    void  setVoxelPosition(const Eigen::Vector3d&);
    void  setInitialPosition(const double &, const double &, const double &);
    void  setInitialPosition(const Eigen::Vector3d &);
    void  setNextDirection(Eigen::Vector3d &);
    void  setRandomInitialPosition(const Eigen::Vector3d &min, const Eigen::Vector3d &max);
    void  setIndex(unsigned int&);

    void setRealPosLog(const Eigen::Vector3d &pos,unsigned t);
    void setRealPosLog(double x, double y, double z, unsigned t);
    void setVoxPosLog(const Eigen::Vector3d &pos,unsigned t);
    void setVoxPosLog(double x, double y, double z, unsigned t);


    void setNumberOfSteps(unsigned T);

    void setColision(unsigned hit_in, unsigned hit_ext, unsigned cross_in, unsigned cross_ext, unsigned t);


};


#endif // WALKER_H
