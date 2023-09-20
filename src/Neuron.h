//!  Class to declare several types of neurons =============================================================/
/*!
*   \details   Class derived from an Obstacle. Defines neurons
*   \author    Inès de Riedmatten
*   \date      February 2023 
*   \version   1.42
=================================================================================================*/


#ifndef NEURON_H
#define NEURON_H

#include "Dendrite.h"

using namespace std;

/// @brief 
class Neuron : public Obstacle
{
public:

    uint8_t nb_dendrites;                   /* Number of dendrites */
    double span_radius;                     /* Radius [mm] inside which all the dendrites are contained */
    std::vector<Dendrite> dendrites;        /* Contains all the dendrites of the neuron*/
    Sphere soma;                    /* Soma of the neuron */
    std::vector<Eigen::Vector2d> Box;       /* Contains the bounding box around the whole neuron */

    /*! Default constructor.*/
    Neuron();
    /*! Constructor
    *  \param dendrites Vector3d<Axon>, dendrites of This.
    *  \param soma Sphere     , soma of This.
    */
    Neuron(std::vector<Dendrite> const& dendrites, Sphere const& soma);
    /*! Constructor
    *  \param dendrites Vector3d<Axon>  , dendrites of This.
    *  \param soma_center Vector3d<Axon>, center of the soma.
    *  \param soma_radius double        , radius of the soma.
    */
    Neuron(std::vector<Dendrite> const& dendrites, Eigen::Vector3d const& soma_center, double const& soma_radius);
    /*! Constructor
    *  \param soma_center Vector3d<Axon>, center of the soma.
    *  \param soma_radius double        , radius of the soma.
    */
    Neuron(Eigen::Vector3d const&soma_center, double const& soma_radius, int const& neuron_id);
    /*! Copy constructor 
    *  \param neuron Neuron
    */
    Neuron(Neuron const& neuron);
    /*! Default destructor.*/
    ~Neuron();
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
    bool checkCollision_branching(Walker &walker, Sphere* const& sphere, Eigen::Vector3d const& step, double const& step_lenght, Collision &colision);
    std::vector<Sphere*> find_neighbor_spheres(Walker &walker, Eigen::Vector3d const& next_step, double const& step_length);
    /**
     * Calculate if position is inside this.
     *
     * @param position Eigen::Vector3d, position of the walker.
     * @param barrier_thickness double, thickness of the cellular barrier.
     * @param swell_              bool, if this swells or not. 
     * @param in_soma_index        int, 0 if in soma, -1 otherwise. 
     * @param in_dendrite_index    int, id of the dendrite it is in, -1 if outside dendrite. 
     * @param in_subbranch_index   int, id of the dendrite subbranch it is in, -1 if outside dendrite. 
     * @return                    bool, true if position is inside this.
     */
    bool isPosInsideNeuron(Eigen::Vector3d const& position,  double const& barrier_thickness, bool const& swell_, int& in_soma_index, 
                           int& in_dendrite_index, int& in_subbranch_index, std::vector<int>& in_sph_index);
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
     * Add dendrite to dendrites. 
     *
     * @param dendrite_to_add Dendrite.
     */
    void add_dendrite(Dendrite& dendrite_to_add);
    // std::vector <int> closest_subbranch(Eigen::Vector3d const& position, int const& dendrite_id, int const& subbranch_id, double const& step_length);
    std::vector <double> get_Volume() const;
    /**
     * Calculates if there is/are intersection(s) between the sphere s and a walker
     * starting at traj_orig, with a direction step_dir. 
     * There can be none, one or two intersections.
     * 
     * Taken from : https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
     *
     * @param intercept1         double, distance from traj_orig to the first potential intercept.
     * @param intercept2         double, distance from traj_orig to the second potential intercept.
     * @param s          Sphere, sphere with which to calculate the intersection.
     * @param step_length        double, step length of the walker
     * @param traj_orig Eigen::Vector3d, trajectory origin of the walker
     * @param c                  double, ||traj.orig - s.center||² - s.radius²
     */
    static bool intersection_sphere_vector(double& intercept1, double& intercept2, Sphere const& s, Eigen::Vector3d const& step_dir, 
                                    double const& step_length, Eigen::Vector3d const& traj_orig, double& c, double& a);
    /**
     * Check if the position is in the close vicinity (bounding boxes) of the soma.
     * @return std::tuple<std::string, int>, {'neuron_part', part_id}
     *                                       neuron_part : "soma", "dendrite" or "none",
     *                                       part_id     : soma_id, dendrite_id or -1.
    */
    bool isNearSoma(Eigen::Vector3d const& position,  double const& barrier_thickness) const;
    /**
     * Check if the position is in the close vicinity (bounding boxes) of a denrite.
     * @return std::tuple<std::string, int>, {'neuron_part', dendrite_id}
     *                                       neuron_part : "dendrite" or "none",
     *                                       dendrite_id : dendrite_id or -1.
    */
    std::vector<int> isNearDendrite(Eigen::Vector3d const& position,  double const& barrier_thickness) const;
    std::vector<int> isNearDendrite(Walker const& walker, Eigen::Vector3d const &step_dir, double const &step_lenght, double const& barrier_thickness) const;
    bool isNearNeuron(Eigen::Vector3d const &position, double const &distance_to_be_inside) const;

    void add_projection();

private:

    static int nb_neurons;                  /* Number of neurons in the simulation*/
    /**
     * Vector of all distances between the position pos and the spheres of the this.
     * 
     * @param pos Eigen::Vector3d, position of the walker.
     * @return std::vector<double>, all the distances. 
     */    
    std::vector<double> Distances_to_Spheres(Eigen::Vector3d const& pos) const;
    /**
     * Vector of all distances between the walker w and the spheres of the this.
     * 
     * @param w Walker.
     * @return std::vector<double>, all the distances. 
     */
    std::vector<double> Distances_to_Spheres(Walker const& w) const;
    /**
     * Assign Dendrites to dendrites. 
     *
     * @param Dendrites std::vector<Dendrite>, set of Dendrites.
     */
    void set_dendrites(std::vector<Dendrite> & Dendrites);
    /**
     * Generate a random span radius in [mm], in the interval [lower_bound, upper_bound].
     *
     * @param lower_bound int, lower bound of the interval, in [mm]. 
     * @param upper_bound int, upper bound of the interval, in [mm].
     */
    void generateSpanRadius(double const& lower_bound=0.2, double const& upper_bound=0.5);
    // /**
    //  * Find the closest dendrite from the soma. Indeed, the walker can go from the 
    //  * soma to the dendrites.
    //  *
    //  * @param position Eigen::Vector3d, initial position of the walker.
    //  * @param step_length       double, length of the step of the walker.
    //  * @returns        dendrite_id int, index of the dendrite. 
    //  */
    // std::vector<int> close_dendrites_from_soma(Eigen::Vector3d const& position, double const& step_length);
};


#endif // NEURON_H
