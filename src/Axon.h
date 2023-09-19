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
    //Projections projections;

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

    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double step_lenght, Collision &colision);
    std::vector<int> checkAxisForCollision(Eigen::Vector3d position, double distance_to_be_inside, int axis);
    bool isPosInsideAxon_(Eigen::Vector3d position, double distance_to_be_inside,std::vector<int> &sph_ids);
    //bool isPosInsideAxon(Eigen::Vector3d &position,  double distance_to_be_inside, double max_radius, std::vector<int> &sph_ids);
    bool intersection_sphere_vector(double &t1, double &t2,Sphere &s, Eigen::Vector3d &step, double &step_length, Eigen::Vector3d &pos, double &c);
    void set_spheres(std::vector<Sphere> spheres_to_add);
    bool isNearAxon(Eigen::Vector3d position, double distance_to_be_inside);
    bool isNearAxon(Walker walker, double distance_to_be_inside);
    //void add_projection(Sphere sphere_to_add);
    bool isWalkerInsideAxon(Walker &walker, double distance_to_be_inside);
    /*! \fn  set_prob_crossings
     *  \brief sets the probability of crossing for all spheres
     */
    void set_prob_crossings(double step_length_pref);   
    
    
};


#endif // AXON_H
