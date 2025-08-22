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
        //projections.clear_projections();
        radius = radius_;
    }
    Axon(Axon const &ax);

    inline bool is_empty(const Box& b);
    inline void extend(Box& b, const Eigen::Vector3d& p);
    void build_axon_grid_spheres(const std::vector<Sphere>& spheres,
                                      double cell_size, double pad);
    inline int neighbor_radius_cells(const HashGrid& G, double query_pad);
    inline bool point_in_inflated_aabb(const Eigen::Vector3d& p,
                                   double d);
    void set_spheres(std::vector<Sphere> &spheres_to_add);
    bool is_point_near_axon(const Eigen::Vector3d& p,
                        double d);
    inline bool segment_aabb_intersect(const Eigen::Vector3d& p0,
                                   const Eigen::Vector3d& p1,
                                   const Box& box,
                                   double& tEnter, double& tExit);
    void gather_candidates_DDA(const Eigen::Vector3d& p0,
                           const Eigen::Vector3d& dir_unit,
                           double L,
                           std::vector<int>& out_ids);
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

    bool ensure_same_compartment_at_hit(const Eigen::Vector3d& p0,
                                           const Eigen::Vector3d& dir_unit,
                                           bool start_inside,
                                           double pad,      // use grid.build_pad
                                           double cell,     // grid.cell
                                           double& t_hit,   // in/out
                                           int max_iter);
    
    int occupancy_at_point(const Eigen::Vector3d& p, double margin, const bool& isintra);
    double signed_distance_to_union(const Eigen::Vector3d& p, double margin);
    void set_prob_crossings(double step_length_pref);
    double minDistance(const Walker& w) const;
    bool isPosInsideAxon(const Eigen::Vector3d& p, double margin);
    
};


#endif // AXON_H
