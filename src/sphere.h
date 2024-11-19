//!  Sphere Obstacle Derived Class =============================================================/
/*!
*   \details   Sphere class derived from an Obstacle. Defines sphere of radius R
*   \author    Remy Gardier
*   \date      January 2021
*   \version   0.0
=================================================================================================*/


#ifndef SPHERE_H
#define SPHERE_H

#include "obstacle.h"


class Sphere : public Obstacle
{
public:

    int id;                 /*!< ID of the sphere       */
    Eigen::Vector3d P;      /*!< Center of the sphere   */
    double radius;          /*!< Radius of the sphere   */
    double volume;
    int object_id;          /*!< ID of the object */
    int object_type;        /*!< Type of the object  (0 : axon, 1 : glial)   */
    int branch_id;          /*!< ID of the branch */

    /*!
     *  \brief Default constructor. Does nothing
     */
    Sphere();

    ~Sphere(); 

    /*!
     *  \param P_ Sphere origin
     *  \param radius_ sphere's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    Sphere(int id_, int object_id_, Eigen::Vector3d P_, double radius_, int object_type_, int branch_id_ = -1, double scale = 1):P(P_*scale), radius(radius_*scale){
        id = id_;
        object_id = object_id_;
        object_type = object_type_;
        branch_id = branch_id_;
        volume = 4./3.*M_PI * (radius_*scale) *  (radius_*scale)  *  (radius_*scale);
        // if sphere is part of axon, no branches
        if (object_type == 0){
            branch_id = -1;
        }
    }

    /*!
     *  \param P_ Sphere origin
     *  \param radius_ sphere's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
   Sphere(Sphere const &sph);

    /*! \fn  checkCollision
     *  \param walker, Walker instance in the simulation.
     *  \param 3d step. Is assumed to be normalized.
     *  \param step_length, length used as the maximum step collision distance.
     *  \param collision, Collision instance to save the collision (if any) details.
     *  \return true only if there was a Collision::hit status. \see Collision.
     *  \brief Basic collision function. Returns the if there was any collision on against the obstacle.
     */
    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    /*! \fn  minDistance
     *  \param walker, Walker instance in the simulation.
     *  \brief Returns the minimum distance from the walker to the sphere. Used to set the reachable
     *  sphere that a given walker can reach.
     */
    double minDistance(Walker &w);

    /*! \fn  minDistance
     *  \param O, position in 3d coordinates
     *  \brief Returns the minimum distance from the position to the sphere. Used to set the reachable
     *  sphere that a given walker can reach.
     */

    double minDistance(Eigen::Vector3d O);

private:

    /*! \fn  handleCollition
     *  \param walker, Walker instance in the simulation.
     *  \param collision, Collision instance to save all the information.
     *  \param step, step vector where to move.
     *  \brief Returns true if it was any analytical collision to the infinite plane
     */
    inline bool handleCollition(Walker& walker, Collision &colision, Eigen::Vector3d& step,double& a,double& b, double& c,double& discr,double& step_length);

};

#endif // SPHERE_H
