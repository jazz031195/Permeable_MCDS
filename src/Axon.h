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
    std::vector<Eigen::Vector2d> Box;
    std::vector <int> proximal_branching; /* Contains the parent axon id as well as the other child id at the proximal branching */
    std::vector <int> distal_branching;   /* Contains the child.ren id at the distal branching */

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon(){};

    ~Axon(){
        // spheres.clear();
        //projections.clear_projections();
    };


    Axon(int id_,  Eigen::Vector3d begin_,Eigen::Vector3d end_ , double radius_){

        id = id_;
        begin = begin_;
        end = end_;
        // spheres.clear();
        //projections.clear_projections();
        radius = radius_;

        // proximal_branching.clear();
        // distal_branching.clear();
    }

    Axon(int const& id_, double const& radius_,  Eigen::Vector3d const& begin_, Eigen::Vector3d const& end_, 
         std::vector <int> proximal_branching_, std::vector <int> distal_branching_):
    Axon(id_, begin_, end_, radius_)
    {
        proximal_branching = proximal_branching_;
        distal_branching   = distal_branching_;
    }

    Axon(Axon const& ax);

    bool checkCollision(Walker &walker, Eigen::Vector3d const& step, double const& step_lenght, Collision &colision);
    std::vector<int> checkAxisForCollision(Eigen::Vector3d position, double distance_to_be_inside, int axis);
    bool isPosInsideAxon_(Eigen::Vector3d position, double distance_to_be_inside,std::vector<int> &sph_ids);
    //bool isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside, double max_radius, std::vector<int> &sph_ids);
    bool intersection_sphere_vector(double &t1, double &t2, Sphere &s, Eigen::Vector3d const& step, double const& step_length, Eigen::Vector3d &pos, double &c, double &a);
    void set_spheres(std::vector<Sphere> const& spheres_to_add);
    bool isNearAxon(Eigen::Vector3d position, double distance_to_be_inside);
    bool isNearAxon(Walker walker, double distance_to_be_inside);
    //void add_projection(Sphere sphere_to_add);
    bool isWalkerInsideAxon(Walker &walker, double distance_to_be_inside);
    bool isPosInsideAxon_(Eigen::Vector3d const&position,  double const& distance_to_be_inside, std::vector<int>& sphere_ids, std::vector<double>& distances);
    /* Calculate the volume of the current axon and update its tortuosity */
    double volumeAxon() const;
    /* Calculate the area of the current axon */
    double areaAxon() const;
       
};


#endif // AXON_H
