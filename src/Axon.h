#ifndef AXON_H
#define AXON_H

#include "sphere.h"
#include "obstacle.h"
#include <vector>

using namespace std;

/// @brief 
class Axon : public Obstacle
{
public:
    int id;
    std::vector<Sphere> spheres;
    double radius;
    Eigen::Vector3d begin;
    Eigen::Vector3d end;
    
    struct Box {
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        double z_min;
        double z_max;
    };

    int nbr_boxes = 10;
    std::vector<Box> boxes;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon(){};

    ~Axon(){
        spheres.clear();
        //projections.clear_projections();
    };


    Axon(int id_,  Eigen::Vector3d begin_,Eigen::Vector3d end_ , double radius_){

        id = id_;
        begin = begin_;
        end = end_;
        spheres.clear();
        //projections.clear_projections();
        radius = radius_;
    }
    Axon(Axon const &ax);

    bool checkCollision(const Walker &walker, Eigen::Vector3d &step, const double &step_length, Collision &collision);
    std::vector<int> checkAxisForCollision(const Eigen::Vector3d &position, const double &distance_to_be_inside, const int &axis);
    bool isPosInsideAxon_(const Eigen::Vector3d &position, const double &distance_to_be_inside);
    //bool isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside, double max_radius, std::vector<int> &sph_ids);
    bool intersection_sphere_vector(double &t1, double &t2, const Sphere &s, const Eigen::Vector3d &step, const Eigen::Vector3d &pos);
    void set_spheres(const std::vector<Sphere> &spheres_to_add);
    bool isNearAxon(const Eigen::Vector3d &position, const double &distance_to_be_inside);
    bool isNearAxon(const Walker &walker, const double &distance_to_be_inside);
    //void add_projection(Sphere sphere_to_add);
    bool isWalkerInsideAxon(const Walker &walker,  const double &distance_to_be_inside);
    /*! \fn  set_prob_crossings
     *  \brief sets the probability of crossing for all spheres
     */
    void find_all_intersections(const Walker &walker, const Eigen::Vector3d &step, const double &distance, std::vector<std::pair<double, size_t>> &dist_and_indices);
    bool FindSphereinAxon(const Eigen::Vector3d &position, const double &distance_to_be_inside, std::vector<int> &sph_ids);
    double minDistance(const Walker &w);
    double minDistance(const Eigen::Vector3d &O);
    bool isInsideBox(const Eigen::Vector3d& position, const double &distance_to_be_inside);
    double distanceToBox(const Eigen::Vector3d& O, const Box &box);
    std::vector<int> findCommonIntegers(const std::vector<int>& vec1, const std::vector<int>& vec2, const std::vector<int>& vec3);
    bool isInsideOneBox(const Eigen::Vector3d& position, const double &distance_to_be_inside, const Box &box);
    double distanceToBoxes(const Eigen::Vector3d& O);
    
};


#endif // AXON_H
