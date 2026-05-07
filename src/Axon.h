#ifndef AXON_H
#define AXON_H

#include "sphere.h"
#include "obstacle.h"
#include <vector>
#include <unordered_map>

using namespace std;

/// @brief 
class Axon : public Obstacle
{
public:
    double radius;
    Eigen::Vector3d begin;
    Eigen::Vector3d end;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Axon();

    ~Axon();

    Axon(int id_,  Eigen::Vector3d begin_,Eigen::Vector3d end_ , double radius_){

        id = id_;
        begin = begin_;
        end = end_;
        spheres.clear();
        radius = radius_;
    }
    Axon(Axon const &ax);

    int getObstacleType() const override;
    
};


#endif // AXON_H
