#include "Dendrite.h"

#include <iostream>

#include <cmath>
#include <typeinfo>

using namespace Eigen;
using namespace std;

int Dendrite::nb_dendrites = 0;

Dendrite::Dendrite()
{
    id = nb_dendrites++;
}

Dendrite::Dendrite(std::vector<Axon> const& subbranches_):Dendrite()
{
    subbranches = subbranches_;
}

Dendrite::Dendrite(Dendrite const& dendrite)
{
    subbranches  = dendrite.subbranches;
    Box          = dendrite.Box;
}

Dendrite::~Dendrite()
{
    nb_dendrites--;
}

bool Dendrite::checkCollision(Walker &walker, Vector3d const&step_dir, double const&step_lenght, Collision &collision)
{
    bool bool_collision = false;
    for(size_t i=0; i < subbranches.size(); i++)
    {
        bool_collision = subbranches[i].checkCollision(walker, step_dir, step_lenght, collision);
        if(bool_collision)
            return true;
    }
    return false;
}

bool Dendrite::isPosInsideDendrite(Vector3d const& position,  double const& barrier_thickness, bool const& swell_)
{
    vector<int> dummy;
    for(size_t s=0; s < subbranches.size(); ++s)
    {
        if(subbranches[s].isPosInsideAxon_(position, barrier_thickness, dummy))
            return true;
    }
    return false;
}

double Dendrite::minDistance(Vector3d const& pos) const
{
    // TO IMPLEMENT
    assert(0);
}

double Dendrite::minDistance(Walker const& walker) const
{
    // TO IMPLEMENT
    assert(0);
}

void Dendrite::add_subbranch(Axon& subbranch)
{
    subbranches.push_back(subbranch);
    subbranch.id = subbranches.size() - 1;
}

void Dendrite::set_dendrite(vector<Axon> const& subbranches_)
{
    subbranches = subbranches_;
    add_projection();
}

int Dendrite::get_nb_subbranches()
{
    return subbranches.size();
}

double Dendrite::volumeDendrite() const
{
    double volume = 0;
    for (size_t i=0; i < subbranches.size(); i++)
        volume += subbranches[i].volumeAxon();
    
    return volume;
}

double Dendrite::areaDendrite() const
{
    double area = 0;
    for (size_t i=0; i < subbranches.size(); i++)
        area += subbranches[i].areaAxon();
    
    return area;
}

void Dendrite::add_projection()
{
    for(int axis=0; axis < 3 ; axis++)
    {
        // Contains the minimum axis projection of the subbranches axon_projections
        double min_axis_projection = 1000;
        // Contains the maximum axis projection of the subbranches axon_projections
        double max_axis_projection = 0;

        for(size_t b=0; b < subbranches.size(); b++)
        {            
            // Calculate the box around each subbranch
            int size         = subbranches[b].spheres.size() - 1;
            double center0   = subbranches[b].spheres[0].center[axis];
            double centerEnd = subbranches[b].spheres[size].center[axis];
            // Take the 1st sphere as it is the largest in subbranch 0 (funnel)
            // TODO [ines] : take the max radius of the subbranch
            double radius    = subbranches[b].spheres[0].radius;
            // Take a margin to account for approximation errors
            double eps       = 2.0 * radius;

            // Find the smallest projection of all subbranches
            if(min_axis_projection > center0 - radius - eps)
                min_axis_projection = center0 - radius - eps;
            if(min_axis_projection > centerEnd - radius - eps)
                min_axis_projection = centerEnd - radius - eps;
            // Find the largest projection of all subbranches
            if(max_axis_projection < center0 + radius + eps)
                max_axis_projection = center0 + radius + eps;
            if(max_axis_projection < centerEnd + radius + eps)
                max_axis_projection = centerEnd + radius + eps;
        }
        // The box around the whole dendrite is between the min & max projections
        Box.push_back({min_axis_projection, max_axis_projection});  
    }    
}