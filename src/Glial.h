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
#include <unordered_map>


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
    inline Box make_empty_box() {
        const double inf = std::numeric_limits<double>::infinity();
        return {+inf, -inf, +inf, -inf, +inf, -inf};
    }

    struct HashGrid {
        double cell = 1.0;
        Eigen::Vector3d origin = Eigen::Vector3d::Zero();
        std::vector<std::pair<int,int>> objs;                 // (branch_id, cell_id)
        std::unordered_map<uint64_t, std::vector<int>> buckets;
        Box big_box;
        double build_pad = 0.0;                                // store pad used at build
        double max_radius_plus_pad = 0.0;                     
    };

    HashGrid grid; /*!< grid for fast access to processes */

    Glial();

    ~Glial();

    Glial(int id_, Sphere soma_)
    {
        id = id_;
        soma = soma_;
        processes = {};

    }

    Glial(Glial const &gl);
    void set_spheres(std::vector<Sphere> &spheres_to_add);
    bool checkCollision(const Walker &walker,  Eigen::Vector3d &step, const double& step_lenght, Collision &collision);
    void set_prob_crossings(double step_length_pref);
    double minDistance(const Walker &w) const;
    void build_glia_grid_processes(const std::vector<std::vector<Sphere>>& processes, double cell_size, double pad);
    void set_up_glialcell(std::vector<Sphere> &spheres_to_add);
    bool is_point_near_glia(const Eigen::Vector3d& p, double d);
    inline bool raySphere(const Eigen::Vector3d& p0, const Eigen::Vector3d& dir_unit, const Eigen::Vector3d& C, double R, double& t_enter, double& t_exit);
    void gather_candidates_DDA(const Eigen::Vector3d& p0, const Eigen::Vector3d& dir_unit, double L, std::vector<int>& out_ids);
    inline bool segment_aabb_intersect(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Box& box, double& tEnter, double& tExit);
    inline void extend(Box& b, const Eigen::Vector3d& p);
    inline bool is_empty(const Box& b);
    int isPosInsideGlialCell(const Eigen::Vector3d& p, double margin, const double& L);
    inline bool point_in_inflated_aabb(const Eigen::Vector3d& p, double d);
    int occupancy_at_point(const Eigen::Vector3d& p,
                              double margin,
                              const bool& isintra, const double& L) const;
    inline int neighbor_radius_cells(const HashGrid& G, double query_pad);


};

#endif // GLIAL_H