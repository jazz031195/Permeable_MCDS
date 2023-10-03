//!  Class to declare a dendrite =============================================================/
/*!
*   \details   Class derived from an Obstacle. Defines a dendrite
*   \author    Inès de Riedmatten
*   \date      April 2023 
*   \version   1.42
=================================================================================================*/

#ifndef DENDRITE_H
#define DENDRITE_H

#include "Axon.h"
#include <iostream>
#include <random>
#include <tuple>
#include <string>

/// @brief 
class Dendrite : public Obstacle
{
public:

std::vector<Axon> subbranches;          /* Contains all the subbranches of the dendrite*/
std::vector<Eigen::Vector2d> Box;       /* Contains the bounding box around the whole dendrite */

/*! Default constructor.*/
Dendrite();
/*! Constructor
*  \param subbranches Vector3d<Axon>, subbranches of This.
*/
Dendrite(std::vector<Axon> const& subbranches);
/*! Copy constructor 
*  \param Dendrite dendrite
*/
Dendrite(Dendrite const& dendrite);
/*! Default destructor.*/
~Dendrite();
/**
 * Check if there is a collision with walker and this.
 *
 * @param walker             Walker.
 * @param step_dir  Eigen::Vector3d, direction of the step.
 * @param step_lenght        double, length of the step. 
 * @param collision       Collision, class to handle collision. 
 * @return                     bool, true if collision with this.
 */
bool checkCollision(Walker &walker, Eigen::Vector3d const&step_dir, double const&step_lenght, Collision &collision);
/**
 * Calculate if position is inside this.
 *
 * @param position Eigen::Vector3d, position of the walker.
 * @param barrier_thickness double, thickness of the cellular barrier.
 * @param swell_              bool, if this swells or not. 
 * @param w                 Walker, if this swells or not. 
 * @return                    bool, true if position is inside this.
 */
bool isPosInsideDendrite(Eigen::Vector3d const& position,  double const& barrier_thickness, bool const& swell_);
/**
 * Minimal distance between a position pos and this.
 *
 * @param pos Eigen::Vector3d, position of the walker.
 * @return minimal distance bewteen `pos` and this.
 */
double minDistance(Eigen::Vector3d const& pos) const;
/**
 * Minimal distance between a walker and this.
 *
 * @param walker .
 * @return minimal distance bewteen `walker` and this.
 */
double minDistance(Walker const& walker) const;
/**
 * Add subbranch to dendrite. 
 *
 * @param subbranch Axon.
 */
void add_subbranch(Axon& subbranch);
/**
 * Get the number of subbranches of a dendrite 
 */
int get_nb_subbranches();
/**
 * Get the volume of a dendrite. 
 */
double volumeDendrite() const;
/**
 * Get the area of a dendrite. 
 */
double areaDendrite() const;
/**
 * Assign subbranches to subbranches. 
 *
 * @param subbranches std::vector<Axon>, set of axons.
 */
void set_dendrite(std::vector<Axon> const& subbranches);
/**
 * Adds the projection (the bounding box) around the whole dendrite.
 * 
 * @param dendrite_id int, id of the dendrite
*/
void add_projection();

private:

    static int nb_dendrites;                  /* Number of dendrites in the simulation*/
    // /**
    //  * Vector of all distances between the position pos and the spheres of the this.
    //  * 
    //  * @param pos Eigen::Vector3d, position of the walker.
    //  * @return std::vector<double>, all the distances. 
    //  */    
    // std::vector<double> Distances_to_Spheres(Eigen::Vector3d const& pos) const;
    // /**
    //  * Vector of all distances between the walker w and the spheres of the this.
    //  * 
    //  * @param w Walker.
    //  * @return std::vector<double>, all the distances. 
    //  */
    // std::vector<double> Distances_to_Spheres(Walker const& w) const;
    // /**
    //  * Check if the position is in the close vicinity (bounding boxes) of this.
    //  * @return std::tuple<std::string, int>, {'neuron_part', part_id}
    //  *                                       neuron_part : "soma", "dendrite" or "none",
    //  *                                       part_id     : soma_id, dendrite_id or -1.
    // */
    // /**
    //  * Calculates if there is/are intersection(s) between the sphere s and a walker
    //  * starting at traj_orig, with a direction step_dir. 
    //  * There can be none, one or two intersections.
    //  * 
    //  * Taken from : https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    //  *
    //  * @param intercept1         double, distance from traj_orig to the first potential intercept.
    //  * @param intercept2         double, distance from traj_orig to the second potential intercept.
    //  * @param s          Dynamic_sphere, sphere with which to calculate the intersection.
    //  * @param step_length        double, step length of the walker
    //  * @param traj_orig Eigen::Vector3d, trajectory origin of the walker
    //  * @param c                  double, ||traj.orig - s.center||² - s.radius²
    //  */
    // bool intersection_sphere_vector(double &intercept1, double &intercept2, Dynamic_Sphere const&s, Eigen::Vector3d const&step_dir, 
    //                                 double const&step_length, Eigen::Vector3d const&traj_orig, double &c) const;
    
};


#endif // DENDRITE_H
