#ifndef BLOOD_VESSEL_H
#define BLOOD_VESSEL_H

#include "obstacle.h"
#include <vector>
#include <unordered_map>

using namespace std;

class Blood_Vessel : public Obstacle
{
    public:
        double pressure_diff; //Pa/mm
        double viscosity; // Pa * s
        double flow;
        double max_velocity; // mm/s
        double min_velocity; // mm/s
        std::vector<Eigen::Vector3d> skeleton;
        double radius;
        bool use_blood_random_direction;
        
        /*! \brief Default constructor. Does nothing
         */
        Blood_Vessel();

        ~Blood_Vessel();

        Blood_Vessel(const Blood_Vessel &bv);

        Blood_Vessel(int id_, double radius_){
        
                id = id_;
                radius = radius_; // mm
                viscosity = 3* 1e-3; // Pa * s for blood at 37 degree C
                pressure_diff = 0.0; 
                flow = 0.0;
                max_velocity = 0.0; // mm/s
                min_velocity = 0.0; // mm/s
        } 

        int getObstacleType() const override;

        void set_spheres(std::vector<Sphere> &spheres_to_add);
                                              
        void distance_to_skeleton(const Walker &w, double& min_dist, Eigen::Vector3d& tangent);
        
        void WalkerVelocity(const Walker &w, double& v, Eigen::Vector3d& flow_direction);

        void set_bv_parameters(const double &pressure_diff_);
        double velocity(const double &radial_distance);

        Eigen::Vector3d generate_random_point_on_sphere(double std);
        Eigen::Vector3d apply_bias_toward_target(const Eigen::Vector3d& point, const Eigen::Vector3d& target);

        Eigen::Matrix3d rotation_matrix_from_vectors(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2);
        Eigen::Vector3d biased_direction_from_tangent(const Eigen::Vector3d& tangent);

        bool checkCollision(const Walker& walker, Eigen::Vector3d& step, const double& step_length, Collision& collision) override;


};

#endif // BLOOD_VESSEL_H