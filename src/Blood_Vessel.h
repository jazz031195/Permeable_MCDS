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
        std::vector<Sphere> spheres; 
        double pressure_diff; //Pa/m
        double viscosity; // Pa * s
        double flow;
        double max_velocity; // mm/s
        std::vector<Eigen::Vector3d> skeleton;
        double radius;

        struct Box {
                double x_min;
                double x_max;
                double y_min;
                double y_max;
                double z_min;
                double z_max;
        };

        inline Box make_empty_box() {
                const double inf = std::numeric_limits<double>::infinity();
                return {+inf, -inf, +inf, -inf, +inf, -inf};
        }

        struct HashGrid {
                double cell = 1.0;
                Eigen::Vector3d origin = Eigen::Vector3d::Zero();
                std::vector<int> objs;                 // (cell_id)
                std::unordered_map<uint64_t, std::vector<int>> buckets;
                Box big_box;
                double build_pad = 0.0;                                // store pad used at build
                double max_radius_plus_pad = 0.0;                     
        };

        HashGrid grid; /*!< grid for fast access to spheres */


        Blood_Vessel();

        ~Blood_Vessel();

        Blood_Vessel(const Blood_Vessel &bv);

        Blood_Vessel(int id_, double radius_, double pressure_diff_){
        
                id = id_;
                radius = radius_;
                viscosity = 3* 1e-3; // Pa * s for blood at 37 degree C
                pressure_diff = pressure_diff_; // Pa/m
                flow = (M_PI*pressure_diff*(radius_*radius_*radius_*radius_))/(8*viscosity); // mm^3/s
                max_velocity = (pressure_diff/(4*viscosity))*(radius*radius); // mm/s
                
        }

        void set_spheres(std::vector<Sphere> &spheres_to_add);
        void build_bv_grid_spheres(const std::vector<Sphere>& spheres_to_add,
                                              double cell_size, double pad);
        inline int neighbor_radius_cells(const HashGrid& G, double query_pad);
        inline bool point_in_inflated_aabb(const Eigen::Vector3d& p,
                                       double d);
        void distance_to_skeleton(const Walker &w, double& min_dist, Eigen::Vector3d& tangent);
        double minDistance(const Walker& w);
        double minDistance(const Eigen::Vector3d& p);
        inline bool is_empty(const Box& b);
        inline void extend(Box& b, const Eigen::Vector3d& p);
        void velocity(const Walker &w, double& v, Eigen::Vector3d& flow_direction);
        inline bool segment_aabb_intersect(const Eigen::Vector3d& p0,
                                   const Eigen::Vector3d& p1,
                                   const Box& box,
                                   double& tEnter, double& tExit);
        void gather_candidates_DDA(const Eigen::Vector3d& p0,
                           const Eigen::Vector3d& dir_unit,
                           double L,
                           std::vector<int>& out_ids);
        bool isPosInsideBlood_Vessel(const Eigen::Vector3d& p, double margin, const double& L);
        int occupancy_at_point(const Eigen::Vector3d& p,
                                  double margin,
                                  const bool& isintra, const double & L) const;
        inline bool raySphere(const Eigen::Vector3d& p0, 
                          const Eigen::Vector3d& dir_unit, // must be unit
                          const Eigen::Vector3d& C,
                          double R,
                          double& t_enter,
                          double& t_exit);
        bool checkCollision(const Walker& walker,
                               Eigen::Vector3d& step,
                               const double& step_length,
                               Collision& collision);



};

#endif // BLOOD_VESSEL_H