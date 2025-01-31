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

    struct Box {
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        double z_min;
        double z_max;
    };

    std::vector<Box> boxes; /*!< boxes around soma and processes */
    Box big_box; /*!< big box around soma and processes */

    Glial();

    ~Glial();

    Glial(int id_, Sphere soma_)
    {
        id = id_;
        soma = soma_;
        processes = {};
        boxes = {};

    }

    Glial(Glial const &gl);
    bool isNearGlialCell(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int> &branches);
    void set_spheres(std::vector<Sphere> &spheres_to_add);
    bool isPosInsideGlialCell(const Eigen::Vector3d& position, const double& distance_to_be_inside);
    std::vector<std::tuple<int, int>> checkAxisForCollision(const Eigen::Vector3d &position, double distance_to_be_inside, int axis, const std::vector<int> &branches);
    bool intersection_sphere_vector(double &t1, double &t2, const Sphere &s, const Eigen::Vector3d &step, const Eigen::Vector3d &pos);
    bool checkCollision(const Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, Collision &collision);
    vector<tuple<int, int>> findCommonIntegers(const vector<vector<tuple<int, int>>>& axisVectors);
    void set_prob_crossings(double step_length_pref);
    void find_all_intersections(const Walker &walker, const Eigen::Vector3d &step, const double &distance, std::vector<std::pair<double, std::tuple<int, int>>> &dist_and_indices);
    bool FindSphereinGlial(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<std::tuple<int, int>> &items);
    double minDistance(const Walker &w);
    bool isInsideBox(const int& i, const Eigen::Vector3d& position, const double &distance_to_be_inside);
    double distanceToBox(const int& i, const Eigen::Vector3d& O);
    double distanceToBigBox(const Eigen::Vector3d& O);
    bool isInsideBigBox(const Eigen::Vector3d &position, const double &distance_to_be_inside);

};

#endif // GLIAL_H