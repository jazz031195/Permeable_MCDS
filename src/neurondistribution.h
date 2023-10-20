//!  NeuronDistribution Class ====================================================================/
/*!
*   \details   Class to construct a substrate of neurons placed in a single voxel structure.
*   \author    Inès de Riedmatten
*   \date      Janvier 2023
*   \version   0.2
=================================================================================================*/

#ifndef NEURONDISTRIBUTION_H
#define NEURONDISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "Neuron.h"


class NeuronDistribution 
{
public:
    std::vector<Neuron> neurons;                    /*!< All neurons of the simulation                                              */
    Eigen::Vector3d min_limits_vx;                  /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits_vx;                  /*!< voxel max limits (if any)                                                  */
    
    struct branching_pt{
        Eigen::Vector3d origin;
        Eigen::Vector3d direction;
        vector<Eigen::Vector3d> children_direction;
        int subbranch_id;
    };

    NeuronDistribution(){}
    /**
     *  @brief Constructor.
     *  @param num_obstacles              int, number of obstacles in the distribution.
     *  @param icvf                    double, intra compartment volume fraction.
     *  @param min_limits_vx_ Eigen::Vector3d, minimum limit of the simulation voxel.
     *  @param max_limits_vx_ Eigen::Vector3d, maximum limit of the simulation voxel.
     *  @param step_length double, length of a walker step, based on Einstein equation.
    */
    NeuronDistribution(int const& num_obstacles, double const& icvf, Eigen::Vector3d const& min_limits_vx_, 
                       Eigen::Vector3d const& max_limits_vx_, double const& step_length);
    /**
     *  @brief Populate the simulation voxel with Neurons.
    */ 
    void createSubstrate();
    /**
     *  @brief Prints the neurons positions in a file or output stream.
     *  @param out std::ostream where to write the info.
    */
    void printSubstrate(std::ostream &out) const;
    void printSubstrate_swc(std::ostream &out) const;
    
private:

    int num_obstacles;                              /*!< number of neurons to fit inside the substrate                              */
    double icvf;                                    /*!< Achieved intra-celular volum fraction in the substrate                     */
    double step_length;                             /*!< Length of a walker step, based on Einstein equation                        */

    struct projection_pt{                           /*!< Structure to calculate the projection of a sphere                          */
        double position;
        int neuron_id;
        int sphere_id;
    };
    std::vector<projection_pt> projections_x;       /*!< Projections of all dendrites spheres on x axis                             */
    std::vector<projection_pt> projections_y;       /*!< Projections of all dendrites spheres on y axis                             */
    std::vector<projection_pt> projections_z;       /*!< Projections of all dendrites spheres on z axis                             */
    /**
     *  Compute the intracompartment volume fraction (ICVF) of the substrate. 
     *  @return ICVF double.
    */
    std::tuple<double, double, double> computeICVF(double const& min_distance_from_border) const;
    /**
     *  Compute the max_limits_vx based on the target icvf and the radiis of the axons 
    */
    // void computeMinimalSize(std::vector<double> const& radiis, double &icvf_, Eigen::Vector3d &l) const;
    bool isSphereColliding(Eigen::Vector3d const& sphere_center, double const& sphere_radius, vector<Eigen::Vector3d> const& soma_centers);
    /**
     * Check if a sphere sph is colliding with this.
     * 
     * @param sph Sphere, sphere to test.
    */
    bool isSphereColliding(Sphere const& sph);
    /**
     * Check if a sphere characterized by its center sphere_center and 
     * its radius sphere_radius is colliding with this.
     * 
     * @param sphere_center Eigen::Vector3d, center of the sphere.
     * @param sphere_radius double         , radius of the sphere.
    */
    bool isSphereColliding(Eigen::Vector3d const& sphere_center, double const& sphere_radius);
    bool isSphereCollidingSphere(Eigen::Vector3d const& pos1, Eigen::Vector3d const& pos2, double const& radius1, double const& radius2, double const& minDistance) const; 
    void createTwinSphere(Eigen::Vector3d &center, double const& sphere_radius, bool &discard_dendrite, size_t const& j);
    /**
     * Check if a position pos is inside the simulation voxel
     * 
     * @param pos       Eigen::Vector3d, position to check.
     * @param distance_to_border double, minimal distance tolerated from the 
     *                                   borders (e.g. to give space for radius).
     * @return bool, true if the position is inside.
    */
    bool isInVoxel(Eigen::Vector3d const& pos, double const& distance_to_border) const;
    /**
     * Grows the dendrites from a neuron.
     * 
     * @param neuron Neuron.
    */
    void growDendrites(Neuron& neuron);
    /**
     * Generate the number of consecutive branching of a dendrite.
     *
     * @param lower_bound int, lower bound of the interval, in [mm]. 
     * @param upper_bound int, upper bound of the interval, in [mm].
     */
    int generateNbBranching(int const& lower_bound=5, int const& upper_bound=7);
    /**
     * Generate the length of a segment/subbranch of the dendrite.
     *
     * @param lower_bound double, lower bound of the interval, in [mm]. 60e-3
     * @param upper_bound double, upper bound of the interval, in [mm]. 100e-3
     */
    double generateLengthSegment(double const& lower_bound=5e-3, double const& upper_bound=15e-3);
    /**
     * Generate the bifurcation angle between two segments/subbranchs at a branching point.
     *
     * @param lower_bound double, lower bound of the interval, in [rad]. 
     * @param upper_bound double, upper bound of the interval, in [rad].
     */
    double generateBifurcationAngle(double const& lower_bound=(M_PI/4 - M_PI/16), double const& upper_bound=(M_PI/4 + M_PI/16));

    Eigen::Vector3d generatePointOnSphere(Eigen::Vector3d const& center, double const& radius) const;
    /**
     * Grow a subbranch/segment of dendrite.
     *
     * @param dendrite                  Dendrite, dendrite to grow. 
     * @param parent  tuple<Eigen::Vector3d,int>, < origin of the subbranch, parent subbranch id>.
     * @param dendrite_direction Eigen::Vector3d, direction of the dendrite / attractor point.
     * @param l_segment                   double, length the segment.
     * @param sphere_radius               double, radius of the spheres in the segment.
     * @param proximal_end          TODO : to complete [ines]
     * @param distal_end
     * @param min_distance_from_border
     */
    branching_pt growSubbranch(Dendrite& dendrite, branching_pt const& parent, double const& l_segment, double const& sphere_radius, 
                               std::vector<int> const& proximal_end, std::vector<int> const& distal_end, double const& min_distance_from_border, 
                               bool& stop_growth, int const& branch_id, Sphere* soma);

    Eigen::Vector3d rotateDirection(Eigen::Vector3d const& direction, double const& angle) const;
    void print_tree(Neuron* const& neuron, std::vector<std::vector<int>> const& nodes) const;
    std::vector<int> find_tree_path(Neuron* const& neuron, int const& dendrite_id, int const& subbranch_id) const;
    std::vector<std::vector<int>> find_tree_paths(Neuron* const& neuron, int const& dendrite_id) const;





};

#endif // NEURONDISTRIBUTION_H
