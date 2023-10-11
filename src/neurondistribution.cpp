#include "neurondistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"

using namespace std;
using namespace Eigen;

NeuronDistribution::NeuronDistribution(int const& num_obstacles_, double const& icvf_, Vector3d const& min_limits_vx_, Vector3d const& max_limits_vx_, double const& step_length_):
min_limits_vx(min_limits_vx_), max_limits_vx(max_limits_vx_), num_obstacles(num_obstacles_), icvf(icvf_), step_length(step_length_)
{
    neurons.clear();
    projections_x.clear();
    projections_y.clear();
    projections_z.clear();
    
    string message = "neurons : " + std::to_string(this->num_obstacles) + " \n";
    SimErrno::info(message, std::cout);
}

// void NeuronDistribution::computeMinimalSize(std::vector<double> const& radiis, double &icvf_, Vector3d &l) const
// {

//     /*A little heuristic for complicated ICVF: > 0.7*/
//     if (icvf_ >= 0.7 && icvf_ < 0.99)
//         icvf_ += 0.01;

//     double area = 0;

//     for (uint i = 0; i < radiis.size(); i++)
//         area += radiis[i] * radiis[i] * M_PI;

//     double l_ = sqrt(area / icvf_);

//     l = {l_, l_, l_};
// }


void NeuronDistribution::createSubstrate()
{

    uint repetition = 1;
    uint max_adjustments = 5;
    // double best_icvf = 0;
    // Eigen::Vector3d best_max_limits;

    bool achieved = false;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);

    uint adjustments = 0;
     // We increease 1% the total area. (Is prefered to fit all the spheres than achieve a perfect ICVF.)
    // double adj_increase = icvf*0.01;
    while(!achieved){

        // double target_icvf = this->icvf+adjustments*adj_increase;
        double soma_radius = 10e-3; //mm
        // Let enough distance for the radius and for a step_length so that 
        // mirroring border conditions are ok
        double min_distance_from_border = barrier_tickness + soma_radius;
        for(uint t = 0 ;  t < repetition; t++){
            neurons.clear();
            vector<Eigen::Vector3d> soma_centers;
            soma_centers.clear();
            
            for(int i = 0 ; i < num_obstacles; i++){
                unsigned stuck = 0;
                while(++stuck <= 10000){
                    // double t = udist(gen);
                    // double x = (t*max_limits_vx[0] + (1-t)*min_limits_vx[0]);
                    // t        = udist(gen);
                    // double y = (t*max_limits_vx[1] + (1-t)*min_limits_vx[1]);
                    // t        = udist(gen);
                    // double z = (t*max_limits_vx[2] + (1-t)*min_limits_vx[2]);
                    double x = (max_limits_vx[0] - min_limits_vx[0])/2;
                    double y = (max_limits_vx[1] - min_limits_vx[1])/2;
                    double z = (max_limits_vx[2] - min_limits_vx[2])/2;

                    Eigen::Vector3d soma_center = {x, y, z};
                   
                    // If too close to the border, discard
                    if ((x - min_limits_vx[0] <= min_distance_from_border) ||
                        (max_limits_vx[0] - x <= min_distance_from_border) ||
                        (y - min_limits_vx[1] <= min_distance_from_border) ||
                        (max_limits_vx[1] - y <= min_distance_from_border) ||
                        (z - min_limits_vx[2] <= min_distance_from_border) ||
                        (max_limits_vx[2] - z <= min_distance_from_border))
                        {
                            continue;
                        } 
                    if (!isSphereColliding(soma_center, soma_radius, soma_centers))
                    {
                        soma_centers.push_back(soma_center);
                        break;
                    }  
                }
            } // end for neurons
            for(size_t i=0; i < soma_centers.size(); ++i)
            {
                Neuron neuron(soma_centers[i], soma_radius, neurons.size());
                growDendrites(neuron); 
                neuron.add_projection();
                neurons.push_back(neuron);
                cout << " End of neuron " << i << endl;
            }


            double icvf, somaFraction, dendritesFraction;
            tie(icvf, somaFraction, dendritesFraction) = computeICVF(0);
            cout << icvf << somaFraction << dendritesFraction << endl;
            achieved = true;

            // if(this->icvf - best_icvf  < 0.0005){
            //     achieved = true;
            //     break;
            // }
        }

        adjustments++;
        if(adjustments > max_adjustments){
            break;
        }
    }

    // spheres = best_spheres;
    // max_limits = best_max_limits;

    // //TODO cambiar a INFO
    // int perc_;
    // double icvf_current = computeICVF(spheres,min_limits_vx, max_limits_vx,perc_);
    // string  message = "Percentage of spheres  selected: "+ to_string(double(perc_)/radiis.size()*100.0)
    //        + "%,\nICVF achieved: " + to_string(icvf_current*100) + "  ("+ to_string( int((icvf_current/icvf*100))) + "% of the desired icvf)\n";
    // SimErrno::info(message,cout);
}

void NeuronDistribution::growDendrites(Neuron& neuron)
{
    // Store all the starting points of dendrites, on the soma of the neuron
    std::vector<Eigen::Vector3d> start_dendrites;
    int max_tries = 10000;

    for(uint8_t i = 0; i < neuron.nb_dendrites; ++i)
    {   
        cout << "dendrite " << i << endl;
        int tries = 0;
        int nb_branching = 3;//generateNbBranching();
        // Radius of each dendrite sphere [mm]
        double sphere_radius = 0.5e-3;
        // Don't initiate dendrite too close from the borders
        double min_distance_from_border = barrier_tickness + sphere_radius + step_length;
        
        while(tries < max_tries)
        {
            Vector3d dendrite_start = generatePointOnSphere(neuron.soma.center, neuron.soma.radius);
            Dendrite dendrite;

            while(!isInVoxel(dendrite_start, min_distance_from_border) || (isSphereColliding(dendrite_start, sphere_radius)))
            {
               dendrite_start = generatePointOnSphere(neuron.soma.center, neuron.soma.radius);
            }

            // If the vector is not already contained in start_dendrites, add it. 
            // Otherwise, decrement i and do one more round
            if((i != 0) && 
               std::count(start_dendrites.begin(), start_dendrites.end(), dendrite_start)){ i--; tries++; }
            else
            {
                start_dendrites.push_back(dendrite_start);
                Eigen::Vector3d dendrite_direction = dendrite_start - neuron.soma.center;
                dendrite_direction.normalize();
                vector<Vector3d> children_dir;
                vector<branching_pt> branching_points     = {{dendrite_start, dendrite_direction, children_dir, 0}};
                vector<branching_pt> branching_points_new;
                int branch_id      = 0;
                int largest_node   = 0;
                // Create the subbranches
                for(int b=0; b < nb_branching; ++b)
                {
                    // Length of a segment before branching
                    double l_segment = 240e-3 / double(nb_branching);//generateLengthSegment();
                    // Number of spheres per segment
                    int nb_spheres   = l_segment / (sphere_radius / 4.0); //Let's assume that dendrites have a radius of 0.5microns so far
                    vector<int> proximal_branch;
                    vector<int> distal_branch;
                    if(b == 0)
                    {
                        proximal_branch = {-1};
                        distal_branch = {largest_node + 1, largest_node + 2};
                        largest_node  = largest_node + 2;
                        bool stop_growth = false;
                        cout << "P id " << branching_points[0].subbranch_id << endl;
                        cout << "C id " << branch_id << endl;
                        cout << "dist " << distal_branch[0] << distal_branch[1] << endl;
                        branching_pt branching_pt_new = growSubbranch(dendrite, branching_points[0], nb_spheres, sphere_radius, proximal_branch, distal_branch, 
                                                                      min_distance_from_border, stop_growth, branch_id, &neuron.soma);
                        branch_id++;
                        branching_points[0] = branching_pt_new;
                        cout << "Subbranch 1 added " << endl;
                        if(stop_growth)
                            break;
                    }
                    else
                    {
                        for(size_t p=0; p < branching_points.size(); p++)
                        {
                                                        // The branching point is split into 2 children
                            for(int c=0; c < static_cast<int>(branching_points[p].children_direction.size()); c++)
                            {
                                proximal_branch = {branching_points[p].subbranch_id, branch_id + 1 - 2*c};
                                if((b + 1) == nb_branching)
                                    distal_branch = {-1};
                                else
                                    distal_branch = {largest_node + 1, largest_node + 2};
                                largest_node  = largest_node + 2;
                                branching_points[p].direction = branching_points[p].children_direction[c];
                                bool stop_growth = false;
                                cout << "P id " << branching_points[p].subbranch_id << endl;
                                cout << "C id " << branch_id << endl;
                                cout << "largest id " << largest_node << endl;
                                cout << "prox " << proximal_branch[0] << proximal_branch[1] << endl;
                                if(distal_branch.size() == 2)
                                    cout << "dist " << distal_branch[0] << distal_branch[1] << endl;
                                branching_pt branching_pt_new = growSubbranch(dendrite, branching_points[p], nb_spheres, sphere_radius, 
                                                                              proximal_branch, distal_branch, min_distance_from_border, 
                                                                              stop_growth, branch_id, &neuron.soma);
                                
                                branch_id++; 
                                if(stop_growth)
                                {
                                    cout << "border reached" << endl;
                                    break;
                                }
                                branching_points_new.push_back(branching_pt_new);  
                            }
                        }

                        branching_points = branching_points_new;
                        branching_points_new = {};
                    }   
                }        
            }
            neuron.add_dendrite(dendrite);
            break;
        }
    } 
}

std::tuple<double, double>  phi_theta_to_target (Eigen::Vector3d parent_dir) 
{
    /*
    Find phi and theta angles in spherical coordinates
    see : http://douglashalse.com/index.php/2019/09/25/spherical-coordinates/
    */
   
    Vector3d vector_to_target = parent_dir.normalized();
    double phi_to_target;
    double theta_to_target;
    // theta is angle between (1,0) and (x,y)
    // varies between 0 and 2pi
    theta_to_target = atan2(vector_to_target[1], vector_to_target[0]);
    if (theta_to_target < 0)
    {
        theta_to_target += 2 * M_PI;
    }
    // phi is angle between (0,0,1) and (x,y,z)
    // varies between 0 and pi
    if (vector_to_target[2] == 0)
    {
        phi_to_target = M_PI / 2.0;
    }
    else if (vector_to_target == Eigen::Vector3d({0, 0, -1}))
    {
        phi_to_target = M_PI;
    }
    else if (vector_to_target == Eigen::Vector3d({0, 0, 1}))
    {
        phi_to_target = 0;
    }
    else
    {
        // varies between -pi/2 and pi/2
        phi_to_target = atan((sqrt(vector_to_target[0] * vector_to_target[0] + vector_to_target[1] * vector_to_target[1])) / vector_to_target[2]);
        if (phi_to_target < 0)
        {
            phi_to_target += M_PI;
        }
    }
    return std::make_tuple(phi_to_target, theta_to_target);
}

NeuronDistribution::branching_pt NeuronDistribution::growSubbranch(Dendrite& dendrite, NeuronDistribution::branching_pt const& parent, 
                                      int const& nb_spheres, double const& sphere_radius, vector<int> const& proximal_end, 
                                      vector<int> const& distal_end, double const& min_distance_from_border, bool& stop_growth,
                                      int const& branch_id, Sphere* soma)
{
    Eigen::Vector3d begin;
    Axon subbranch(branch_id, sphere_radius, begin, begin, proximal_end, distal_end);
    std::vector<Sphere> spheres_to_add;
    spheres_to_add.clear();

    Eigen::Vector3d center = {0, 0, 0};
    bool discard_dendrite  = false;
   
    for(int sphere_id=0; sphere_id < nb_spheres ; ++sphere_id)
    {
        if(branch_id == 0)
            center = double(sphere_id) * parent.direction * sphere_radius / 4.0 + parent.origin;
        else
            center = double(sphere_id + 1) * parent.direction * sphere_radius / 4.0 + parent.origin;

        if(isInVoxel(center, min_distance_from_border))
        {
            if (!isSphereColliding(center, sphere_radius))
            {
                Sphere sphere_to_add(sphere_id, branch_id, center, sphere_radius);
                if(sphere_id == 0)
                {
                    if(proximal_end[0] == -1)
                    {
                        sphere_to_add.add_neighbor(soma);
                        soma->add_neighbor(new Sphere(sphere_to_add));
                    }  
                    else
                    {
                        sphere_to_add.add_neighbor(new Sphere(dendrite.subbranches[proximal_end[0]].spheres[dendrite.subbranches[proximal_end[0]].spheres.size()-1]));
                        dendrite.subbranches[proximal_end[0]].spheres[dendrite.subbranches[proximal_end[0]].spheres.size()-1].add_neighbor(new Sphere(sphere_to_add));
                        cout << "b id " << branch_id << endl;
                        cout << "current sph id" << sphere_to_add.id << " ax id " << sphere_to_add.ax_id << endl;
                        cout << "sph id " << dendrite.subbranches[proximal_end[0]].spheres[dendrite.subbranches[proximal_end[0]].spheres.size()-1].id << endl;
                        cout << "ax id " << dendrite.subbranches[proximal_end[0]].spheres[dendrite.subbranches[proximal_end[0]].spheres.size()-1].ax_id<< endl;
                        if(branch_id > proximal_end[1])
                        {
                            sphere_to_add.add_neighbor(new Sphere(dendrite.subbranches[proximal_end[1]].spheres[0]));
                            dendrite.subbranches[proximal_end[1]].spheres[0].add_neighbor(new Sphere(sphere_to_add));
                            cout << "sph id " << dendrite.subbranches[proximal_end[1]].spheres[0].id<< endl;
                            cout << "ax id " << dendrite.subbranches[proximal_end[1]].spheres[0].ax_id<< endl;
                        }
                    }
                }
                else
                {
                    spheres_to_add[spheres_to_add.size() - 1].add_neighbor(new Sphere(sphere_to_add));
                    sphere_to_add.add_neighbor(new Sphere(spheres_to_add[spheres_to_add.size() - 1]));
                }
                
                spheres_to_add.push_back(sphere_to_add);
            }
        }
        else
        {
            if(spheres_to_add.size() > 0)
            {
                subbranch.set_spheres(spheres_to_add);
                dendrite.add_subbranch(subbranch);
            }
            return {};
        } 
    }
    if (!discard_dendrite)
    {
        subbranch.set_spheres(spheres_to_add);
        dendrite.add_subbranch(subbranch);
    } 
    
    float phi_to_target, theta_to_target;
    tie(phi_to_target, theta_to_target) = phi_theta_to_target(parent.direction);

    random_device dev;
    mt19937 rng(dev());
    normal_distribution<float> theta_distr(M_PI / 4.0, M_PI / 16.0);
    float delta_theta = theta_distr(rng) / 2.0;
    float theta = theta_to_target + delta_theta;
    float delta_x, delta_y, delta_z;
    vector<Vector3d> children_dir;
    for(int c=0; c < 2; c++)
    {
        delta_x = cos(theta);
        delta_y = sin(theta);
        delta_z = parent.direction[2];

        Vector3d new_dir(delta_x, delta_y, delta_z);
        new_dir = new_dir.normalized();
        children_dir.push_back(new_dir);

        // phi = phi;
        theta = theta_to_target - delta_theta;
    }
    
    uniform_real_distribution<double> rot_distr(0, M_PI);

    double rot = rot_distr(rng);
    double ux = parent.direction[0];
    double uy = parent.direction[1];
    double uz = parent.direction[2];
    Matrix3d m;
    m << cos(rot) + ux*ux*(1 - cos(rot)), ux*uy*(1-cos(rot)) - uz*sin(rot), ux*uz*(1-cos(rot)) + uy*sin(rot),
        uy*ux*(1-cos(rot)) + uz*sin(rot), cos(rot) + uy*uy*(1 - cos(rot)), uy*uz*(1-cos(rot)) - ux*sin(rot),
        uz*ux*(1-cos(rot)) - uy*sin(rot), uz*uy*(1-cos(rot)) + ux*sin(rot), cos(rot) + uz*uz*(1 - cos(rot)); 
    

    children_dir[0] = (m * children_dir[0]).normalized();
    children_dir[1] = (m * children_dir[1]).normalized();
    // cout << acos(children_dir[0].dot(children_dir[1])) << endl;
    // cout << 2*delta_theta << endl;
    // cout << subbranch.projections.axon_projections[0][0] << " " << subbranch.projections.axon_projections[0][1];
    // cout << subbranch.projections.axon_projections[1][0] << " " << subbranch.projections.axon_projections[1][1];
    // cout << subbranch.projections.axon_projections[2][0] << " " << subbranch.projections.axon_projections[2][1];
    
    // Return the next branching point
    return {spheres_to_add[spheres_to_add.size()-1].center, parent.direction, children_dir, branch_id};
    
}



double NeuronDistribution::generateLengthSegment(double const& lower_bound, double const& upper_bound)
{
    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<> segmentLength(lower_bound, upper_bound);

    return segmentLength(rng);
}


double NeuronDistribution::generateBifurcationAngle(double const& lower_bound, double const& upper_bound)
{
    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<double> bifurcationAngle(lower_bound, upper_bound);

    return bifurcationAngle(rng);
}

Vector3d NeuronDistribution::generatePointOnSphere(Vector3d const& center, double const& radius) const
{
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::normal_distribution<double> distribution(0.0, 1.0);

    // Find a point at the surface of the soma
    double x = distribution(generator); 
    double y = distribution(generator);
    double z = distribution(generator);

    // Avoid division by 0 while normalizing
    while(x==0 && y==0 && z==0)
    {
        x = distribution(generator);
        y = distribution(generator);
        z = distribution(generator);
    }
    double normalization_factor = sqrt(x*x + y*y + z*z);
    x = x/normalization_factor*radius + center[0];
    y = y/normalization_factor*radius + center[1];
    z = z/normalization_factor*radius + center[2];
    
    return {x, y, z};   
}

Vector3d NeuronDistribution::rotateDirection(Vector3d const& direction, double const& angle) const
{

    // Compute a unit vector v1 that makes the desired angle with "direction"
    // Define a random unit vector u
    Vector3d u;
    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<double> generator(-1, 1);
    for (int i=0; i < u.size(); i++)
        u[i] = generator(rng);
    u /= sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);

    // v1 has an angle "angle" with direction
    Vector3d v1;
    for (int i=0; i < u.size(); i++)
        v1[i] = cos(angle)*direction[i] + sin(angle)*u[i];
    v1 /= sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    
    return v1;
}

void NeuronDistribution::createTwinSphere(Vector3d &center, double const& sphere_radius, bool &discard_dendrite, size_t const& j)
{
    Vector3d new_center = center;
    for (size_t axis=0; axis < 3; ++axis)
    {
        double distance_from_min_lim = center[axis] - min_limits_vx[axis];
        if (distance_from_min_lim < sphere_radius)
        {
            new_center[axis] = max_limits_vx[axis] + distance_from_min_lim;
        }
        double distance_from_max_lim = max_limits_vx[axis] - center[axis];
        if (distance_from_max_lim < sphere_radius)
        {
            new_center[axis] = min_limits_vx[axis] - distance_from_max_lim;
        }
    }
    if(!isSphereColliding(new_center, sphere_radius))
    {
        discard_dendrite = false;
        center = new_center;
    }
}

bool NeuronDistribution::isSphereColliding(Sphere const& sph) 
{
    Vector3d position = sph.center;
    double distance_to_be_inside = sph.radius + 2 * barrier_tickness;
    int dummy;
    vector<int> dummy2;
    for (unsigned i = 0; i < neurons.size() ; i++){
        bool isinside = neurons[i].isPosInsideNeuron(position, distance_to_be_inside, dummy, dummy, dummy, dummy2);
        if (isinside)
            return true;
    }
    return false;
}

bool NeuronDistribution::isSphereColliding(Vector3d const& sphere_center, double const& sphere_radius, vector<Vector3d> const& soma_centers) 
{
    double distance_to_be_inside = 2 * barrier_tickness;
    for (unsigned i = 0; i < soma_centers.size() ; i++){
        bool isinside = isSphereCollidingSphere(sphere_center, soma_centers[i], sphere_radius, sphere_radius, distance_to_be_inside);
        if (isinside)
            return true;
    }
    return false;
}

bool NeuronDistribution::isSphereCollidingSphere(Vector3d const& pos1, Vector3d const& pos2, double const& radius1, double const& radius2, double const& minDistance) const 
{
    Vector3d m = pos1 - pos2;
    double distance_to_sphere = m.norm() - radius1 - radius2;

    return distance_to_sphere < minDistance;
}

bool NeuronDistribution::isSphereColliding(Vector3d const& sphere_center, double const& sphere_radius) 
{
    double distance_to_be_inside = sphere_radius + 2 * barrier_tickness;
    int dummy;
    vector<int> dummy2;
    for (unsigned i = 0; i < neurons.size() ; i++){
        bool isinside = neurons[i].isPosInsideNeuron(sphere_center, distance_to_be_inside, dummy, dummy, dummy, dummy2);
        if (isinside)
            return true;
    }
    return false;
}

void NeuronDistribution::printSubstrate(ostream &out) const
{
    out << 1 << endl; //scale
    out << 0 << endl; //volume_inc_perc
    out << 0 << endl; //dyn_perc
    out << icvf << endl;
    out << min_limits_vx[0] << endl; //min_limits [mm]
    out << max_limits_vx[0] << endl; //max_limits [mm]

    for (unsigned i = 0; i < neurons.size(); i++)
    {
        // Print for soma : x y z r bool_active
        // bool_active = 1 if the sphere can be activated (swollen)
        out << neurons[i].soma.center[0] << " " 
        << neurons[i].soma.center[1] << " "
        << neurons[i].soma.center[2] << " "
        << neurons[i].soma.radius    << endl; 
        out << "Soma " + to_string(i) << endl;

        for (size_t j = 0; j < neurons[i].dendrites.size(); j++)
        {
            // vector<double> min = {10, 10, 10};
            // vector<double> max = {0, 0, 0};
            for (size_t k = 0; k < neurons[i].dendrites[j].subbranches.size(); k++)
            {
                for (size_t l = 0; l < neurons[i].dendrites[j].subbranches[k].spheres.size(); l++)
                {
                    // Print for each dendrite, each sphere
                    out << neurons[i].dendrites[j].subbranches[k].spheres[l].center[0] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[1] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[2] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].radius << endl; 
                    // for(int axis=0; axis < 3; axis++)
                    // {
                    //     if(neurons[i].dendrites[j].subbranches[k].spheres[l].center[axis] < min[axis])
                    //         min[axis] = neurons[i].dendrites[j].subbranches[k].spheres[l].center[axis];
                    //     if(neurons[i].dendrites[j].subbranches[k].spheres[l].center[axis] > max[axis])
                    //         max[axis] = neurons[i].dendrites[j].subbranches[k].spheres[l].center[axis];
                    // }
                        
                }
                //TODO [ines] : add proximal & distal branching writing
                out << "Segment " + to_string(k) << endl;
                
            }
            out << "Dendrite " + to_string(j) << endl;
            out << neurons[i].dendrites[j].Box[0][0] << " "
                    << neurons[i].dendrites[j].Box[0][1] << " "
                    << neurons[i].dendrites[j].Box[1][0] << " "
                    << neurons[i].dendrites[j].Box[1][1] << " "
                    << neurons[i].dendrites[j].Box[2][0] << " "
                    << neurons[i].dendrites[j].Box[2][1] << endl;
            // cout << "min " << min[0] << " " << min[1] << " " << min[2] << endl;
            // cout << "max " << max[0] << " " << max[1] << " " << max[2] << endl;
        }
        out << "Neuron " + to_string(i) << endl;
    }
}

void NeuronDistribution::printSubstrate_swc(ostream &out) const
{
    /* 
    SWC format :
    SampleID TypeID x y z r ParentID

    SampleID : Sample identifier. A positive integer.
    TypeID : Type identifier. The basic set of types used in NeuroMorpho.org SWC files is:

    -1  - root
     0  - undefined
     1  - soma
     2  - axon
     3  - (basal) dendrite
     4  - apical dendrite
     5+ - custom

    In addition, some SWC-variants use the following types 5 and 6:

    5 - branch point (redundant: branch is a point with multiple children)
    6 - end point (redundant: end point is a point with zero children)

    x : X-position in micrometers
    y : Y-position in micrometers
    z : Z-position in micrometers
    r : Radius in micrometers (half the cylinder thickness)
    ParentID : Parent sample identifier. This defines how points are connected to 
               each other. In a tree, multiple points can have the same ParentID. 
               The first point in the file must have a ParentID equal to -1, which 
               represents the root point. Parent samples must be defined before they 
               are being referred to. By counting how many points refer to the a 
               given parent, the number of its children can be computed. 

    Ref : https://neuroinformatics.nl/swcPlus/ 
    */
    for (size_t i = 0; i < neurons.size(); i++)
    {
        int line  = 1;
    
        // SampleID TypeID x y z r ParentID
        cout << line                 << " " 
        << 1                        << " "
        << neurons[i].soma.center[0] << " "
        << neurons[i].soma.center[1] << " "
        << neurons[i].soma.center[2] << " "
        << neurons[i].soma.radius    << " "
        << -1                       << endl; 
        line++;
        int dendrite_id  = 0;
        vector<vector<int>> path = find_tree_paths(new Neuron(neurons[i]), dendrite_id);
        for(size_t i=0; i < path.size(); ++i)
            for(size_t j=0; j < path[i].size(); ++j)
                cout << path[i][j] << endl;
        // for(size_t i=0; i < (*neuron).dendrites.size(); ++i)
        // {
        //     auto sphere_node = (*neuron).dendrites[i].subbranches[0].spheres[(*neuron).dendrites[i].subbranches[0].spheres.size() - 1];
        //      // SampleID TypeID x y z r ParentID
        //     out << 1                    << " " 
        //     << 3                        << " "
        //     << sphere_node.center[0]    << " "
        //     << sphere_node.center[1]    << " "
        //     << sphere_node.center[2]    << " "
        //     << sphere_node.radius       << " "
        //     << line                     << endl; 
        //     line++;
        //     for(size_t j=0; j < (*neuron).dendrites[i].subbranches[0].distal_branches.size(); ++j)
        //     {
        //         auto next_subbranch = (*neuron).dendrites[i].subbranches[0].distal_branches[]
        //     } 
        // }
    }
}


vector<vector<int>> NeuronDistribution::find_tree_paths(Neuron* const& neuron, int const& dendrite_id) const
{
    int id_subbranch = (*neuron).dendrites[dendrite_id].subbranches.size() - 1;
    vector<vector<int>> paths;
    vector<int> path;
    while(true)
    {
        if((*neuron).dendrites[dendrite_id].subbranches[id_subbranch].distal_branching[0] == -1)
        {
            paths.push_back(find_tree_path(neuron, dendrite_id, id_subbranch));
            id_subbranch--;
        }
        else
            return paths;
    }
}

vector<int> NeuronDistribution::find_tree_path(Neuron* const& neuron, int const& dendrite_id, int const& subbranch_id) const
{
    vector<int> path;
    path.insert(path.begin(), subbranch_id);
    int subbranch_id_tmp = subbranch_id;
    while(true)
    {
        vector<int> proximal_branching = (*neuron).dendrites[dendrite_id].subbranches[subbranch_id_tmp].proximal_branching;
        if(proximal_branching.size() == 2)
        {
            path.insert(path.begin(), proximal_branching[0]);
            subbranch_id_tmp = proximal_branching[0];
        }
        else
            return path;
    }
}

bool NeuronDistribution::isInVoxel(Eigen::Vector3d const& pos, double const& distance_to_border) const
{
    Eigen::Vector3d new_min_limits_vx = {min_limits_vx[0] + distance_to_border, min_limits_vx[1] + distance_to_border, min_limits_vx[2] + distance_to_border};
    Eigen::Vector3d new_max_limits_vx = {max_limits_vx[0] - distance_to_border, max_limits_vx[1] - distance_to_border, max_limits_vx[2] - distance_to_border};

    if ((pos[0] - new_min_limits_vx[0]) < 0 || 
        (pos[1] - new_min_limits_vx[1]) < 0 || 
        (pos[2] - new_min_limits_vx[2]) < 0)
        return false;
    else if ((pos[0] - new_max_limits_vx[0]) > 0 || 
             (pos[1] - new_max_limits_vx[1]) > 0 || 
             (pos[2] - new_max_limits_vx[2]) > 0) 
        return false;

    return true;   
}

tuple<double, double, double> NeuronDistribution::computeICVF(double const& min_distance_from_border) const
{

    if (neurons.size() == 0)
        return make_tuple(0, 0, 0);

    double VolumeVoxel = (max_limits_vx[0] - min_limits_vx[0] - min_distance_from_border) * (max_limits_vx[1] - min_limits_vx[1] - min_distance_from_border) * (max_limits_vx[2] - min_limits_vx[2] - min_distance_from_border);
    double VolumeSoma = 0;
    double VolumeDendrites = 0;
    vector<double> volumeNeuron;
   
    for (size_t i = 0; i < neurons.size(); i++)
    {
        volumeNeuron     = neurons[i].get_Volume();
        VolumeSoma      += volumeNeuron[0];
        VolumeDendrites += volumeNeuron[1];
    }
    double somaFraction      = VolumeSoma / VolumeVoxel;
    double dendritesFraction = VolumeDendrites/ VolumeVoxel;
    double ICVF              = somaFraction + dendritesFraction;
    return make_tuple(ICVF, somaFraction, dendritesFraction);
}
