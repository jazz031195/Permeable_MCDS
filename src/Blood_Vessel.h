#ifndef BLOOD_VESSEL_H
#define BLOOD_VESSEL_H

#include "sphere.h"
#include "cylinder.h"
#include "obstacle.h"
#include <vector>
#include <unordered_map>

using namespace std;

class Blood_Vessel : public Obstacle
{
    public:
        int id;
        Cylinder vessel;
        double pressure_diff; //Pa
        double viscosity; // Pa * s
        Eigen::Vector3d direction_flow; 
        double flow;
        double max_velocity; // mm/s


        Blood_Vessel();

        ~Blood_Vessel();

        Blood_Vessel(const Blood_Vessel &bv);


        Blood_Vessel(int id_,  Eigen::Vector3d begin_,Eigen::Vector3d end_ , double radius_, double flow_){
        
                id = id_;
                vessel = Cylinder(id, begin_, end_, radius_);
                direction_flow = (end_ - begin_).normalized();
                viscosity = 3* 1e-3; // Pa * s for blood at 37 degree C
                flow = flow_;
                pressure_diff = (8*viscosity*flow)/(M_PI*vessel.radius*vessel.radius*vessel.radius*vessel.radius); // Hagen–Poiseuille equation
                max_velocity = (pressure_diff/(4*viscosity))*(vessel.radius*vessel.radius); // mm/s
        }

        double velocity(const Walker &w);
        double minDistance(const Eigen::Vector3d& O);
        bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);
        bool isPosInsideBloodVessel(const Eigen::Vector3d& p, double margin);
        double minDistance(const Walker &w);

};

#endif // BLOOD_VESSEL_H