//!  Glial Obstacle Derived Class =============================================================/
/*!
*   \details   Glial class derived from an Obstacle
*              in the direction set by begin, end.
*   \author    Jasmine Nguyen-Duc
*   \date      February 2024
*   \version   1.42
=================================================================================================*/

#ifndef GLIAL_H
#define GLIAL_H

#include "sphere.h"
#include "obstacle.h"
#include <vector>


using namespace std;

/// @brief
class Glial : public Obstacle
{
    public : 
    int id;                                         /*!< ID of glial */
    Sphere soma;                                    /*!< soma of glial */
    std::vector<std::vector<Sphere>> processes; /*!< ramification spheres of glial */
    std::vector<Eigen::Vector2d> Box;               /*!< Box with <min, max> for each axis (x,y,z) */
    std::vector<std::vector<Eigen::Vector2d>> Box_branch; /*!< Box with <min, max> for each axis (x,y,z) */

    Glial();

    ~Glial();

    Glial(int id_, Sphere soma_)
    {
        id = id_;
        soma = soma_;
        // intialise box to soma
        // create box around that one sphere
        init_box();

    }

    Glial(Glial const &gl);
    void init_box();
    bool isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside);
    bool isNearDendrite(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int>& dendrite_ids);
    void set_spheres(const std::vector<Sphere> &spheres_to_add);
    void set_box_neuron();
    bool isPosInsideGlialCell(const Eigen::Vector3d& position, const double& distance_to_be_inside);
    bool isPosInsideGlialCell_verbose(const Eigen::Vector3d& position, const double& distance_to_be_inside);
    std::vector<int> checkAxisForCollision_soma(Eigen::Vector3d position, double distance_to_be_inside, int axis);
    std::vector<int> checkAxisForCollision_dendrite(Eigen::Vector3d position, double distance_to_be_inside, int axis, int dendrite_to_check);
    bool intersection_sphere_vector(double &t1, double &t2, Sphere &s, Eigen::Vector3d &step, const double &step_length, const Eigen::Vector3d &pos);
    bool checkCollision(Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, Collision &colision);
    std::vector<int> findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3);
    void set_prob_crossings(double step_length_pref);
    bool isPosInsideGlialCell_(const Eigen::Vector3d& position, const double& distance_to_be_inside);
    void find_all_intersections(const Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, std::vector<double>& dist_intersections, std::vector<Sphere*>& spheres);
    bool FindSphereinGlial(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<Sphere*> &sph);
    double minDistance(Walker &w);

};

#endif // GLIAL_H