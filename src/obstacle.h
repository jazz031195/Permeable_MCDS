//!  Obstacle Base Class ==============================================================================/
/*!
*   \details   Father class to define the base of any other obstacle (wall or substrate)
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42
 =====================================================================================================*/

#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "collision.h"
#include "sphere.h"
#include "walker.h"
#include <Eigen/Core>
#include <random>
#include <unordered_map>

using namespace std;

class Obstacle
{
public:
    int id;                         /*!< Unique id of the simulation                                                */
    std::vector<Sphere> spheres; 
    int count_perc_crossings;       /*!< Auxiliar value to count the number of percolatin crossings in a simulation */
    double percolation;             /*!< Percolation value between 0 and 1.                                         */
    double prob_cross_e_i;
    double prob_cross_i_e;
    double diffusivity_i;
    double diffusivity_e;
    

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

    HashGrid grid;
    

    /*! \fn  Obstacle
     *  \brief Default constructor. Does nothing.
     */


    Obstacle();

    Obstacle(const Obstacle& obs);

    virtual int getObstacleType() const = 0;

    void elasticBounceAgainsPlane(Eigen::Vector3d& ray_origin, Eigen::Vector3d& normal, double& t, Eigen::Vector3d &step);

    double clamp(double value, double lower, double upper);

    inline bool is_empty(const Box& b);

    inline void extend(Box& b, const Eigen::Vector3d& p);

    void set_spheres(std::vector<Sphere> &spheres_to_add);

    void build_bv_grid_spheres(const std::vector<Sphere>& spheres_to_add, double cell_size, double pad);

    inline bool raySphere(const Eigen::Vector3d& p0, const Eigen::Vector3d& dir_unit, const Eigen::Vector3d& C, double R, double& t_enter, double& t_exit);

    bool isPosInsideObstacle(const Eigen::Vector3d& p, double margin);

    void gather_candidates_AABB(const Eigen::Vector3d& p0, const Eigen::Vector3d& dir, double L, const HashGrid& grid, std::vector<int>& out_ids);

    bool checkCollision(const Walker& walker, Eigen::Vector3d& step, const double& step_length, Collision& collision);

    double minDistance(const Walker& w);

    double minDistance(const Eigen::Vector3d& p);

    void setPercolation(double& percolation_);

    void setDiffusion(double& diffusivity_i_, double& diffusivity_e_);

    void setProbabilities(double &prob_cross_e_i_, double &prob_cross_i_e_);

    void set_prob_crossings(double step_length_pref);


};

#endif // OBSTACLE_H
