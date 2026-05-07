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

    Glial();

    ~Glial();

    Glial(int id_)
    {
        id = id_;

    }

    Glial(Glial const &gl);

    int getObstacleType() const override;
};

#endif // GLIAL_H