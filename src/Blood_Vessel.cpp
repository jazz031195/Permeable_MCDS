
#include "Blood_Vessel.h"
#include "Eigen/Dense"
#include <Eigen/Geometry>
#include <Eigen/Core>
#include "constants.h"
#include <numeric>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <set>
#include <limits>
#include <cmath>
#include <limits>
#include <algorithm> // for std::min initializer_list

using namespace Eigen;
using namespace std;

Blood_Vessel::Blood_Vessel()
{}

Blood_Vessel::~Blood_Vessel()
{}

Blood_Vessel::Blood_Vessel(const Blood_Vessel &bv)
{
    id = bv.id;
    vessel = bv.vessel;
    pressure_diff = bv.pressure_diff;
    viscosity = bv.viscosity;
    direction_flow = bv.direction_flow;
    flow = bv.flow;
    max_velocity = bv.max_velocity;
};

double Blood_Vessel::velocity(const Walker &w){

    Eigen::Vector3d pos = w.pos_v;

    double distance = minDistance(pos);

    if (distance >= 0) {
        return 0.0;
    }

    double radial_distance = abs(vessel.radius + distance); // distance is negative inside the vessel

    double v = (pressure_diff/(4*viscosity))*(vessel.radius*vessel.radius - radial_distance*radial_distance); // Parabolic profile

    return v; // mm/s

}

double Blood_Vessel::minDistance(const Eigen::Vector3d& O){

    return vessel.minDistance(O);

}

double Blood_Vessel::minDistance(const Walker &w){

    Eigen::Vector3d O;
    w.getVoxelPosition(O);
    return vessel.minDistance(O);

}

bool Blood_Vessel::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision){
    
    return vessel.checkCollision(walker, step, step_lenght, colision);
}

bool Blood_Vessel::isPosInsideBloodVessel(const Eigen::Vector3d& p, double margin){

    double distance = minDistance(p);

    if (distance < margin) {
        return true;
    }
    else{
        return false;
    }

}