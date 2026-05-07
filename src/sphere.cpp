#include "sphere.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>
#include <random>

using namespace Eigen;
using namespace std;

Sphere::Sphere()
{}

Sphere::~Sphere()
{}

Sphere::Sphere(const Sphere &sph)
{

    P               = sph.P;
    radius          = sph.radius;
    id              = sph.id;
    volume          = sph.volume;
    object_id       = sph.object_id;
    object_type     = sph.object_type;
    branch_id       = sph.branch_id;
    percolation = sph.percolation;
    prob_cross_e_i = sph.prob_cross_e_i;
    prob_cross_i_e = sph.prob_cross_i_e;
    diffusivity_i = sph.diffusivity_i;
    diffusivity_e = sph.diffusivity_e;
    count_perc_crossings = sph.count_perc_crossings;

}

bool Sphere::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{
    //Origin of the ray
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d m = O - P;

    //distance to the sphere center.
    double distance_to_sphere = m.norm();
    double d_ = distance_to_sphere - radius;

    //If the minimum distance from the walker to the sphere is more than
    // the actual step size, we can discard this collision.
    if(d_> EPS_VAL){
        if(d_ > step_lenght+barrier_tickness){
            return false;
        }
    }

    double a = 1 ; 
    double b = step.dot(m) ; 
    double c = distance_to_sphere*distance_to_sphere - radius*radius;

    double discr = b*b - a*c;

    //No real roots
    if(discr < 0.0){
        colision.type = Collision::null;
        return false;
    }

    //if we arrived here we need to compute the quadratic equation.
    return handleCollition(walker,colision,step,a,b,c,discr,step_lenght);

}

inline bool Sphere::handleCollition(Walker& walker, Collision &colision, Vector3d& step,double& a,double& b, double& c,double& discr,double& step_length){

    double t1 = (-b - sqrt(discr))/a;

    double t2 = (-b + sqrt(discr))/a;


    //if we are completely sure that no collision happened
    if( ( (t1 < 0.0) || (t1 > step_length+barrier_tickness) ) && ( (t2 < 0.0) || (t2 > step_length+barrier_tickness)) ){
        colision.type = Collision::null;
        return false;
    }

    // a spin that's bouncing ignores collision at 0 (is in a wall)
    if(walker.status == Walker::bouncing){

        //if the collision are too close or negative.
        if( ( (t1 < EPS_VAL) || (t1 > step_length)) && (( t2 < EPS_VAL) || (t2 > step_length)) ){
            colision.type = Collision::null;
            return false;
        }

        if( t1 >= EPS_VAL && t1 < t2)
            colision.t = fmin(t1,step_length);
        else
            colision.t = fmin(t2,step_length);
    }
    else{
        if( t1>0.0 && t1 <t2)
            colision.t = fmin(t1,step_length);
        else
            colision.t = fmin(t2,step_length);
    }

    colision.type = Collision::hit;
    colision.obstacle_ind = id;

    if(c<-1e-10){
        colision.col_location = Collision::inside;
        walker.in_obj_index = id;
    }
    else if(c>1e-10){
        colision.col_location = Collision::outside;
    }
    else{
        colision.col_location = Collision::unknown;
    }

    colision.rn = c;

    colision.collision_point = walker.pos_v + colision.t*step;
    colision.obstacle_type = sph_obstacle_type;
    

    // Membrane permeability    
    if((this->percolation>0.0)){
        if(colision.type == Collision::hit && colision.col_location != Collision::voxel){

            std::mt19937 gen_perm;          // Random engine for permeability
            std::random_device rd;
            gen_perm.seed(rd());
            std::uniform_real_distribution<double> udist(0,1);
            
            double _percolation_ = udist(gen_perm); 

            double dynamic_percolation = 0.0;
            
            if (colision.col_location == Collision::inside){ 
                dynamic_percolation =  this->prob_cross_i_e; 
            } 

            else if (colision.col_location == Collision::outside){
                dynamic_percolation = this->prob_cross_e_i;
            } 

            if( dynamic_percolation - _percolation_ > EPS_VAL ){            
                count_perc_crossings++;
                colision.perm_crossing      = _percolation_;
                colision.bounced_direction  = step; 
                return true;
            }
        }  
    }

    colision.perm_crossing = 0.;


    if (fabs(a) < EPS_VAL){
        colision.col_location = Collision::on_edge;
        colision.bounced_direction = -step;
    }
    else{

        /* For a sphere, normal direction is equal to colision point */
        //Normal point
        Eigen::Vector3d normal = (colision.collision_point - P).normalized();

        Eigen::Vector3d temp_step = step;
        elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);

        colision.bounced_direction = temp_step.normalized();

    }
    
    return true;

}

double Sphere::minDistance(Walker &w){

    //Origin of the ray
    Vector3d O;
    w.getVoxelPosition(O);
    Vector3d m = O - P;
    // minimum distance to the sphere center.
    double distance_to_sphere = m.norm();

    //Minimum distance to the sphere wall.
    double d_ = (distance_to_sphere - radius);
    return d_;

}

double Sphere::minDistance(Eigen::Vector3d O){

    Vector3d m = O - P;
    // minimum distance to the sphere center.
    double distance_to_sphere = m.norm();

    //Minimum distance to the sphere wall.
    double d_ = (distance_to_sphere - radius);
    return d_;

}

void Sphere::setPercolation(double &percolation_)
{
    percolation = percolation_;
}

void Sphere::setDiffusion(double &diffusivity_i_, double &diffusivity_e_){
    diffusivity_i = diffusivity_i_;
    diffusivity_e = diffusivity_e_;
}


void Sphere::setProbabilities(double &prob_cross_e_i_, double &prob_cross_i_e_)
{
    prob_cross_e_i = prob_cross_e_i_;
    prob_cross_i_e = prob_cross_i_e_;
}

void Sphere::set_prob_crossings(double step_length_pref){

    double prob_cross_i_e_, prob_cross_e_i_;
    double dse, dsi;

    if (percolation > 0.0){
        
        dse = sqrt(step_length_pref*this->diffusivity_e);
        dsi = sqrt(step_length_pref*this->diffusivity_i);

        prob_cross_i_e_ = percolation * dsi * 2. / 3. / this->diffusivity_i;
        prob_cross_e_i_ = percolation * dse * 2. / 3. / this->diffusivity_e; 

        this->prob_cross_e_i = prob_cross_e_i_ / (1.+ 0.5 * (prob_cross_e_i_ + prob_cross_i_e_));
        this->prob_cross_i_e = prob_cross_i_e_ / (1.+ 0.5 * (prob_cross_e_i_ + prob_cross_i_e_));
    }
            

}


void Sphere::elasticBounceAgainsPlane(Eigen::Vector3d &ray_origin, Eigen::Vector3d &normal, double &t, Eigen::Vector3d &step)
{

    Eigen::Vector3d ray =  (-t*step).normalized();//
    double rn = ray.dot(normal);

    //std::cout << "rn in elastic plane :" << rn << std::endl;

    // Caso 3) ni cerca ni paralela
    step = -ray + 2.0*normal*rn;

    //step = (rn>0.0)?normal:(-normal);

}