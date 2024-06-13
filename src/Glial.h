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
    std::vector<Sphere> processes; /*!< ramification spheres of glial */
    std::vector<Eigen::Vector2d> Box;               /*!< Box with <min, max> for each axis (x,y,z) */

    Glial();

    ~Glial();

    Glial(int id_, Sphere soma_)
    {
        id = id_;
        soma = soma_;
        // intialise box to soma
        // create box around that one sphere
        double sph_highest_x_val = soma.P[0]+ soma.radius;
        double sph_lowest_x_val = soma.P[0] -soma.radius;
        double sph_highest_y_val = soma.P[1] +soma.radius;
        double sph_lowest_y_val = soma.P[1] -soma.radius;
        double sph_highest_z_val = soma.P[2] +soma.radius;
        double sph_lowest_z_val = soma.P[2] -soma.radius;
        //x
        Box.push_back({sph_lowest_x_val, sph_highest_x_val});
        //y
        Box.push_back({sph_lowest_y_val, sph_highest_y_val});
        //z
        Box.push_back({sph_lowest_z_val, sph_highest_z_val});

    }

    Glial(Glial const &gl);
    bool isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside);
    void set_spheres(const std::vector<Sphere> &spheres_to_add);
    bool isPosInsideGlialCell(const Eigen::Vector3d& position, const double& distance_to_be_inside);
    std::vector<int> checkAxisForCollision(Eigen::Vector3d position, double distance_to_be_inside, int axis);
    bool intersection_sphere_vector(double &t1, double &t2, Sphere &s, Eigen::Vector3d &step, const double &step_length, const Eigen::Vector3d &pos);
    bool checkCollision(Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, Collision &colision);
    std::vector<int> findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3);
    void set_prob_crossings(double step_length_pref);
    bool isPosInsideGlialCell_(const Eigen::Vector3d& position, const double& distance_to_be_inside);
    void find_all_intersections(const Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, std::vector<double>& dist_intersections, std::vector<int>& spheres_ids);

};

#endif // GLIAL_H