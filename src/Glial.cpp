#include "Glial.h"
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

Glial::Glial()
{}

Glial::~Glial()
{}

Glial::Glial(const Glial &gl)
{
    id = gl.id;
    spheres = gl.spheres;
    grid = gl.grid;
    percolation = gl.percolation;
    prob_cross_e_i = gl.prob_cross_e_i;
    prob_cross_i_e = gl.prob_cross_i_e;
    diffusivity_i = gl.diffusivity_i;
    diffusivity_e = gl.diffusivity_e;
    count_perc_crossings = gl.count_perc_crossings;

};

int Glial::getObstacleType() const { 
    return glial_obstacle_type; 
}


