#include "Axon.h"
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

Axon::Axon()
{}

Axon::~Axon()
{}

Axon::Axon(const Axon &ax)
{
    id = ax.id;
    spheres = ax.spheres;
    begin = ax.begin;
    end = ax.end;
    grid = ax.grid;
    percolation = ax.percolation;
    prob_cross_e_i = ax.prob_cross_e_i;
    prob_cross_i_e = ax.prob_cross_i_e;
    diffusivity_i = ax.diffusivity_i;
    diffusivity_e = ax.diffusivity_e;
    count_perc_crossings = ax.count_perc_crossings;

};

int Axon::getObstacleType() const { 
    return axon_obstacle_type; 
}

