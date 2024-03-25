#include "doctest.h"

#include "neurondistribution.h"
#include <algorithm> // std::sort
#include <random>
#include "simerrno.h"

using namespace std;
using namespace Eigen;    

NeuronDistribution::NeuronDistribution(int const& num_obstacles_, double const& icvf_, Vector3d const& min_limits_vx_, Vector3d const& max_limits_vx_, int const& sphere_overlap_, double const& step_length_):
min_limits_vx(min_limits_vx_), max_limits_vx(max_limits_vx_), num_obstacles(num_obstacles_), sphere_overlap(sphere_overlap_), icvf(icvf_), step_length(step_length_)
{
    neurons.clear();
    projections_x.clear();
    projections_y.clear();
    projections_z.clear();
    
    string message = "neurons : " + std::to_string(this->num_obstacles) + " \n";
    SimErrno::info(message, std::cout);
}

TEST_CASE("updateLUT_add")
{
    int num_obstacles = 1;
    double ICVF = 0;
    Vector3d min_limits_vx = {0.0, 0.0, 0.0};
    Vector3d max_limits_vx = {0.1, 0.1, 0.1};
    int sphere_overlap = 4;
    double step_length = 0.5e-3;
    NeuronDistribution neuronDist(num_obstacles, ICVF, min_limits_vx, max_limits_vx, sphere_overlap, step_length);
    
    SUBCASE("Soma")
    {
        Vector3d sphere_center = {0.05, 0.05, 0.05};
        double sphere_radius   = 10e-3;

        neuronDist.updateLUT(sphere_center, sphere_radius);

        // TODO [ines] : check if starts at left, middle, right part of the subvoxel, check if size correct...

        CHECK(!neuronDist.LUT[0][0][0]);
        CHECK(neuronDist.LUT[81][100][100]);
        CHECK(!neuronDist.LUT[80][80][100]);
    }
    SUBCASE("Soma-dendrites-y")
    {
        Vector3d soma_center = {0.05, 0.05, 0.05};
        double soma_radius     = 10e-3;

        neuronDist.updateLUT(soma_center, soma_radius);

        double dendrite_radius = 0.5e-3;
        Vector3d dendrite_center;
        Vector3d dir = {0.0, 1.0, 0.0};
        for(size_t i=0; i < 50; ++i)
        {
            dendrite_center = soma_center + dir * (soma_radius + double(i) * dendrite_radius / 2.0);
            neuronDist.updateLUT(dendrite_center, dendrite_radius);
        }

        CHECK(!neuronDist.LUT[0][0][0]);
        CHECK(neuronDist.LUT[100][120][100]);
        CHECK(!neuronDist.LUT[80][80][100]);
    }
    SUBCASE("Soma-dendrites-yz")
    {
        Vector3d soma_center = {0.05, 0.05, 0.05};
        double soma_radius     = 10e-3;

        neuronDist.updateLUT(soma_center, soma_radius);

        double dendrite_radius = 0.5e-3;
        Vector3d dendrite_center;
        Vector3d dir = {0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0)};
        for(size_t i=0; i < 50; ++i)
        {
            dendrite_center = soma_center + dir * (soma_radius + double(i) * dendrite_radius / 2.0);
            neuronDist.updateLUT(dendrite_center, dendrite_radius);
        }
        // TODO [ines] : check if starts at left, middle, right part of the subvoxel, check if size correct...
        // imshow((unsigned char *)&neuronDist.LUT[100], SIZE_LUT, SIZE_LUT, 1);
        // show(); 

        CHECK(!neuronDist.LUT[0][0][0]);
        CHECK(neuronDist.LUT[100][120][120]);
        CHECK(!neuronDist.LUT[80][80][100]);
    }
}

void NeuronDistribution::updateLUT(Vector3d const& sphere_center, double const& sphere_radius)
{
    // TODO [ines] : bitwise operations => 8x more resolution
    // uint8_t array[10];
    // constexpr uint8_t size = sizeof(array[0]);
    // size_t index = 75;
    // bool current_state = array[index / size] & (1 << (index % size)) > 0;
    uint16_t min_x, min_y, min_z, max_x, max_y, max_z;
    min_x = floor((sphere_center[0] - sphere_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    assert(floor((sphere_center[0] - sphere_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT)) >= 0);
    max_x = floor((sphere_center[0] + sphere_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    min_y = floor((sphere_center[1] - sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    assert(floor((sphere_center[1] - sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT)) >= 0);
    max_y = floor((sphere_center[1] + sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    min_z = floor((sphere_center[2] - sphere_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
    assert(floor((sphere_center[1] - sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT)) >= 0);
    max_z = floor((sphere_center[2] + sphere_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));

    // cout << min_x << " " << max_x << " " << min_y << " " << max_y << " " << min_z << " " << max_z << endl;

    for(size_t x=min_x; x <= max_x; ++x)
    {
        for(size_t y=min_y; y <= max_y; ++y)
        {
            for(size_t z=min_z; z <= max_z; ++z)
            {
                double x_real = x * (max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT;
                double y_real = y * (max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT;
                double z_real = z * (max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT;
                if(((sphere_center[0] - x_real)*(sphere_center[0] - x_real) + 
                   (sphere_center[1] - y_real)*(sphere_center[1] - y_real) +
                   (sphere_center[2] - z_real)*(sphere_center[2] - z_real)) <= sphere_radius * sphere_radius)
                   {
                        LUT[x][y][z] = true; 
                        // cout << x << " " << y << " " << z << " " << endl;
                   } 
            }
        }
    }
    // imshow((unsigned char *)&LUT[int((min_x + max_x)/2)], SIZE_LUT, SIZE_LUT, 1);
    // show(); 
}

TEST_CASE("updateLUT_remove")
{
    int num_obstacles = 1;
    double ICVF = 0;
    Vector3d min_limits_vx = {0.0, 0.0, 0.0};
    Vector3d max_limits_vx = {0.1, 0.1, 0.1};
    int sphere_overlap = 4;
    double step_length = 0.5e-3;
    NeuronDistribution neuronDist(num_obstacles, ICVF, min_limits_vx, max_limits_vx, sphere_overlap, step_length);
    
    SUBCASE("Soma-dendrites-y")
    {
        vector<Sphere> spheres;
        Vector3d soma_center = {0.05, 0.05, 0.05};
        double soma_radius   = 10e-3;
        Sphere soma(0, 0, soma_center, soma_radius);
        spheres.push_back(soma);

        neuronDist.updateLUT(soma_center, soma_radius);

        double dendrite_radius = 0.5e-3;
        Vector3d dendrite_center;
        Vector3d dir = {0.0, 1.0, 0.0};
        for(size_t i=0; i < 50; ++i)
        {
            dendrite_center = soma_center + dir * (soma_radius + double(i) * dendrite_radius / 2.0);
            Sphere s(0, 0, dendrite_center, dendrite_radius);
            spheres.push_back(s);
            neuronDist.updateLUT(dendrite_center, dendrite_radius);
        }

        CHECK(!neuronDist.LUT[0][0][0]);
        CHECK(neuronDist.LUT[100][144][100]);
        CHECK(!neuronDist.LUT[80][80][100]);

        // Remove 25 spheres
        neuronDist.updateLUT(spheres, 25);

        CHECK(!neuronDist.LUT[0][0][0]);
        CHECK(!neuronDist.LUT[100][144][100]);
        CHECK(!neuronDist.LUT[80][80][100]);
    }
}

void NeuronDistribution::updateLUT(vector<Sphere>& spheres_to_add, int const& nb_spheres_to_remove)
{
    // Remove the last spheres
    for(size_t sph_id = spheres_to_add.size() - 1; sph_id >= spheres_to_add.size() - nb_spheres_to_remove; sph_id--)
    {
        uint16_t min_x, min_y, min_z, max_x, max_y, max_z;
        min_x = floor((spheres_to_add[sph_id].center[0] - spheres_to_add[sph_id].radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
        max_x = floor((spheres_to_add[sph_id].center[0] + spheres_to_add[sph_id].radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
        min_y = floor((spheres_to_add[sph_id].center[1] - spheres_to_add[sph_id].radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
        max_y = floor((spheres_to_add[sph_id].center[1] + spheres_to_add[sph_id].radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
        min_z = floor((spheres_to_add[sph_id].center[2] - spheres_to_add[sph_id].radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
        max_z = floor((spheres_to_add[sph_id].center[2] + spheres_to_add[sph_id].radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));

        for(size_t x=min_x; x < max_x; ++x)
        {
            for(size_t y=min_y; y < max_y; ++y)
            {
                for(size_t z=min_z; z < max_z; ++z)
                {
                    double x_real = x * (max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT;
                    double y_real = y * (max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT;
                    double z_real = z * (max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT;
                    if(((spheres_to_add[sph_id].center[0] - x_real)*(spheres_to_add[sph_id].center[0] - x_real) + 
                    (spheres_to_add[sph_id].center[1] - y_real)*(spheres_to_add[sph_id].center[1] - y_real) +
                    (spheres_to_add[sph_id].center[2] - z_real)*(spheres_to_add[sph_id].center[2] - z_real)) <= spheres_to_add[sph_id].radius * spheres_to_add[sph_id].radius)
                    {
                            LUT[x][y][z] = false;         
                    } 
                }
            }
        }
    }
    
    // Add again the last sphere to keep
    size_t sph_id = spheres_to_add.size() - 1 - nb_spheres_to_remove;
    uint16_t min_x, min_y, min_z, max_x, max_y, max_z;
    min_x = floor((spheres_to_add[sph_id].center[0] - spheres_to_add[sph_id].radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    max_x = floor((spheres_to_add[sph_id].center[0] + spheres_to_add[sph_id].radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    min_y = floor((spheres_to_add[sph_id].center[1] - spheres_to_add[sph_id].radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    max_y = floor((spheres_to_add[sph_id].center[1] + spheres_to_add[sph_id].radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    min_z = floor((spheres_to_add[sph_id].center[2] - spheres_to_add[sph_id].radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
    max_z = floor((spheres_to_add[sph_id].center[2] + spheres_to_add[sph_id].radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));

    for(size_t x=min_x; x < max_x; ++x)
    {
        for(size_t y=min_y; y < max_y; ++y)
        {
            for(size_t z=min_z; z < max_z; ++z)
            {
                double x_real = x * (max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT;
                double y_real = y * (max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT;
                double z_real = z * (max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT;
                if(((spheres_to_add[sph_id].center[0] - x_real)*(spheres_to_add[sph_id].center[0] - x_real) + 
                (spheres_to_add[sph_id].center[1] - y_real)*(spheres_to_add[sph_id].center[1] - y_real) +
                (spheres_to_add[sph_id].center[2] - z_real)*(spheres_to_add[sph_id].center[2] - z_real)) <= spheres_to_add[sph_id].radius * spheres_to_add[sph_id].radius)
                {
                        LUT[x][y][z] = true;  
                        
                } 
            }
        }
    } 

    for(size_t i=0; i < nb_spheres_to_remove; i++)
    {
        spheres_to_add.pop_back();
    }
}

void NeuronDistribution::createSubstrate()
{
    uint repetition    = 1;
    double soma_radius = 10e-3; //mm
    // Let enough distance for the radius and for a step_length so that 
    // mirroring border conditions are ok
    double min_distance_from_border = 10 * barrier_tickness + soma_radius;

    bool achieved = false;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> udistx(min_limits_vx[0] + min_distance_from_border, max_limits_vx[0] - min_distance_from_border);
    uniform_real_distribution<double> udisty(min_limits_vx[1] + min_distance_from_border, max_limits_vx[1] - min_distance_from_border);
    uniform_real_distribution<double> udistz(min_limits_vx[1] + min_distance_from_border, max_limits_vx[1] - min_distance_from_border);

    uint adjustments = 0;
    while(!achieved){

        
        for(uint t = 0 ;  t < repetition; t++)
        {
            neurons.clear();
            vector<Vector3d> soma_centers;
            vector<double> soma_radii;
            soma_centers.clear();

            for(int i = 0 ; i < num_obstacles; i++)
            {
                unsigned stuck = 0;
                while(++stuck <= 10000)
                {
                    double x = udistx(gen);
                    double y = udisty(gen);
                    double z = udistz(gen);

                    Vector3d soma_center = {x, y, z};

                    if (!isSphereColliding(soma_center, soma_radius))
                    {
                        cout << "soma " << i << endl;
                        soma_centers.push_back(soma_center);
                        soma_radii.push_back(soma_radius);
                        updateLUT(soma_center, soma_radius);
                        break;
                    }  
                }
            } // end for neurons
            for(int i = 0 ; i < 200; i++)
            {
                soma_radius = 0.6e-3;
                unsigned stuck = 0;
                while(++stuck <= 10000)
                {
                    double x = udistx(gen);
                    double y = udisty(gen);
                    double z = udistz(gen);

                    Vector3d soma_center = {x, y, z};

                    if (!isSphereColliding(soma_center, soma_radius))
                    {
                        cout << "soma " << i << endl;
                        soma_centers.push_back(soma_center);
                        soma_radii.push_back(soma_radius);
                        updateLUT(soma_center, soma_radius);
                        break;
                    }  
                }
            } // end for neurons
            for(size_t i=0; i < soma_centers.size(); ++i)
            {
                Neuron neuron(soma_centers[i], soma_radii[i], neurons.size());
                growDendrites(neuron, soma_centers, soma_radii); 
                // neuron.add_projection();
                neurons.push_back(neuron);
                cout << " End of neuron " << i << endl;
            }


            double icvf, somaFraction, dendritesFraction;
            tie(icvf, somaFraction, dendritesFraction) = computeICVF(0);
            cout << "ICVF : " << icvf << "soma fraction : " << somaFraction << "dendrites fraction : " << dendritesFraction << endl;
            achieved = true;
        }
    }
}

void NeuronDistribution::growDendrites(Neuron& neuron, vector<Vector3d> soma_centers, vector<double> const& soma_radii)
{
    // Store all the starting points of dendrites, on the soma of the neuron
    std::vector<Eigen::Vector3d> start_dendrites;
    int max_tries = 1000;

    for(uint8_t i = 0; i < neuron.nb_dendrites; ++i)
    {   
        cout << "dendrite " << i << endl;
        int tries = 0;
        int nb_branching = 1;//generateNbBranching();
        // Radius of each dendrite sphere [mm]
        double sphere_radius = 0.6e-3;
        // Don't initiate dendrite too close from the borders
        double min_distance_from_border = barrier_tickness + sphere_radius;
        
        while(tries < max_tries)
        {
            Vector3d dendrite_dir     = generatePointOnSphere();
            Vector3d dendrite_start   = dendrite_dir * neuron.soma.radius + neuron.soma.center;
            Vector3d dendrite_to_test = dendrite_dir * 50e-3 + neuron.soma.center;
            Dendrite dendrite;

            while((!isInVoxel(dendrite_to_test, min_distance_from_border) || 
                  isSphereColliding(dendrite_start, sphere_radius, neuron.soma.center, neuron.soma.radius) ||
                  !avoid_somas(50e-3, soma_centers, soma_radii, neuron.soma.center, dendrite_dir, dendrite_start)) &&
                  (tries < max_tries))
            {
               dendrite_dir     = generatePointOnSphere();
               dendrite_to_test = dendrite_dir * 50e-3 + neuron.soma.center;
               dendrite_start   = dendrite_dir * neuron.soma.radius + neuron.soma.center;
               tries++;
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
                vector<branching_pt> branching_points     = {{neuron.soma.center, dendrite_direction, children_dir, 0, neuron.soma.radius}};
                vector<branching_pt> branching_points_new;
                int branch_id      = 0;
                int largest_node   = 0;
                // Create the subbranches
                for(int b=0; b < nb_branching; ++b)
                {
                    // Length of a segment before branching
                    double l_segment = 100e-3;//generateLengthSegment();
                    vector<int> proximal_branch;
                    vector<int> distal_branch;
                    if(b == 0)
                    {
                        proximal_branch = {-1};
                        if (nb_branching > 1)
                            distal_branch = {largest_node + 1, largest_node + 2};
                        else
                            distal_branch = {-1};
                        largest_node  = largest_node + 2;
                        bool stop_growth = false;
                        // cout << "P id " << branching_points[0].subbranch_id << endl;
                        // cout << "C id " << branch_id << endl;
                        // if (nb_branching > 1)
                        //     cout << "dist " << distal_branch[0] << distal_branch[1] << endl;
                        branching_pt branching_pt_new = growSubbranch(dendrite, branching_points[0], l_segment, sphere_radius, proximal_branch, distal_branch, 
                                                                      min_distance_from_border, branch_id);
                        branch_id++;
                        branching_points[0] = branching_pt_new;

                    } // branching 0
                    else
                    {   
                        cout << "s2" << branching_points.size() << endl;

                        for(size_t p=0; p < branching_points.size(); p++)
                        {
                            // The branching point is split into 2 children
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
                                branching_pt branching_pt_new = growSubbranch(dendrite, branching_points[p], l_segment, sphere_radius, 
                                                                              proximal_branch, distal_branch, min_distance_from_border, 
                                                                              branch_id);
                                
                                branch_id++; 
                                branching_points_new.push_back(branching_pt_new);  
                            }
                        } // grow each branch points

                        branching_points = branching_points_new;
                        branching_points_new = {};
                    } // Subbranchings
                }        
            }
            if(dendrite.subbranches.size() == 0)
                tries++;
            else if(dendrite.subbranches[0].spheres.size() > 50)
            {
                neuron.add_dendrite(dendrite);
                break;
            }
            else
                tries++;
            
        }
    } 
}

tuple<double, double> phi_theta_to_target(Vector3d parent_dir) 
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
        theta_to_target += 2 * M_PI;
    // phi is angle between (0,0,1) and (x,y,z)
    // varies between 0 and pi
    if (vector_to_target[2] == 0)
        phi_to_target = M_PI / 2.0;
    else if (vector_to_target == Eigen::Vector3d({0, 0, -1}))
        phi_to_target = M_PI;
    else if (vector_to_target == Eigen::Vector3d({0, 0, 1}))
        phi_to_target = 0;
    else
    {
        // varies between -pi/2 and pi/2
        phi_to_target = atan((sqrt(vector_to_target[0] * vector_to_target[0] + vector_to_target[1] * vector_to_target[1])) / vector_to_target[2]);
        if (phi_to_target < 0)
            phi_to_target += M_PI;
    }
    return make_tuple(phi_to_target, theta_to_target);
}

NeuronDistribution::branching_pt NeuronDistribution::growSubbranch(Dendrite& dendrite, NeuronDistribution::branching_pt const& parent, 
                                      double const& l_segment, double const& sphere_radius, vector<int> const& proximal_end, 
                                      vector<int> const& distal_end, double const& min_distance_from_border,
                                      int const& branch_id)
{
    Vector3d begin;
    Axon subbranch(branch_id, sphere_radius, begin, begin, proximal_end, distal_end);
    vector<Sphere> spheres_to_add;
    spheres_to_add.clear();

    Vector3d center, center_to_test;
    double last_radius;
    Vector3d end  = parent.origin + (l_segment + parent.radius) * parent.direction;
    int sphere_id = 0;
    float phi_to_target, theta_to_target;
    tie(phi_to_target, theta_to_target) = phi_theta_to_target(parent.direction);
    Vector3d last_sphere_center = parent.origin + parent.radius * parent.direction;
    // size_t nb_spheres_to_interpolate = sphere_overlap / 2 - 1;
    int nb_trials = 0;
    while((nb_trials < 100) & ((last_sphere_center - end).norm() > 1e-6))
    {
        Vector3d next_direction = generateNextDirection(phi_to_target, theta_to_target, 0.03);
        // 1st sphere of 1st subbranch
        if((branch_id == 0) && (sphere_id == 0))
        {
            center         = last_sphere_center;
            center_to_test = parent.origin;
            last_radius    = 10e-3;

        }
        else
        {
            center = next_direction * sphere_radius / 2.0 + last_sphere_center;
            if(sphere_id == 0)
                center_to_test = parent.origin;
            else
                center_to_test = last_sphere_center;

            last_radius = sphere_radius;
        }

        if(isInVoxel(center, min_distance_from_border))
        {
            if (!isSphereColliding(center, sphere_radius, center_to_test, last_radius))
            {
                Sphere sphere_to_add(sphere_id, branch_id, center, sphere_radius);
                if(spheres_to_add.size() > 1)
                    if((sphere_to_add.center - spheres_to_add[spheres_to_add.size() - 1].center).norm() > 0.0004)
                    {
                        cout << "norm " << (sphere_to_add.center - spheres_to_add[spheres_to_add.size() - 1].center).norm() << endl;
                        cout << (center_to_test - center).norm() << endl;
                        cout << "center " << center << endl;
                        cout << "CTT " << center_to_test << endl;
                        // cout << "last " << last_sphere_center << endl;
                        cout << "s -1 " << spheres_to_add[spheres_to_add.size() - 1].center << endl;
                        cout << "s " << sphere_to_add.center << endl;
                    }
                // cout << "p 2 " << sphere_to_add.center << endl;
                spheres_to_add.push_back(sphere_to_add);
                // cout << "p 2 " << spheres_to_add.back().center << endl;
                // if(spheres_to_add.size() > 1)
                //     cout << "p 2 " << (spheres_to_add[spheres_to_add.size() - 2].center - spheres_to_add[spheres_to_add.size() - 1].center).norm() << endl;


                sphere_id++; 
                last_sphere_center = spheres_to_add.back().center;
                updateLUT(center, sphere_radius);
            }
            // Try other directions
            else
            {
                nb_trials++;
                constexpr int nb_spheres_to_remove = 5;
                if(spheres_to_add.size() > nb_spheres_to_remove)
                    updateLUT(spheres_to_add, nb_spheres_to_remove);
                else if(spheres_to_add.size() >= 2)
                    updateLUT(spheres_to_add, spheres_to_add.size() - 1);

                if(spheres_to_add.size() > 1)
                    last_sphere_center = spheres_to_add.back().center;
                else
                    last_sphere_center = parent.origin + parent.radius * parent.direction;
                bool achieved = false;
                int number_dir = 10;
                while (number_dir > 0)
                {
                    next_direction = generateNextDirection(phi_to_target, theta_to_target, 0.03);
                    if(spheres_to_add.size() > 0)
                        center = next_direction * spheres_to_add[spheres_to_add.size() - 1].radius / 2.0 + last_sphere_center; 
                    else
                        center = next_direction * last_radius / 2.0 + last_sphere_center;

                    if(isInVoxel(center, min_distance_from_border) && (!isSphereColliding(center, sphere_radius, center_to_test, last_radius)))
                    {
                        Sphere sphere_to_add(sphere_id, branch_id, center, sphere_radius);

                        // if((sphere_to_add.center - spheres_to_add[spheres_to_add.size() - 1].center).norm() > 0.0004)
                        // {
                        //     cout << "norm " << (sphere_to_add.center - spheres_to_add[spheres_to_add.size() - 1].center).norm() << endl;
                        //     cout << (center_to_test - center).norm() << endl;
                        //     cout << "CTT " << center_to_test << endl;
                        //     cout << "last " << last_sphere_center << endl;
                        //     cout << "s -2 " << spheres_to_add[spheres_to_add.size() - 2].center << endl;
                        //     cout << "s -1 " << spheres_to_add.back().center << endl;
                        //     cout << "s " << sphere_to_add.center << endl;
                        // }
                        
                        // cout << "p 1 " << sphere_to_add.center << endl;
                        spheres_to_add.push_back(sphere_to_add);
                        // cout << "p 1 " << spheres_to_add.back().center << endl;
                        // if(spheres_to_add.size() > 1)
                        //     cout << "p 1 " << (spheres_to_add[spheres_to_add.size() - 2].center - spheres_to_add[spheres_to_add.size() - 1].center).norm() << endl;



                        sphere_id++;
                        achieved = true;
                        last_sphere_center = spheres_to_add.back().center;
                        updateLUT(center, sphere_radius);
                        break;
                    }
                    number_dir--;
                }
                
                if(!achieved)
                {
                    double sphere_radius_shrink = sphere_radius;
                    while (sphere_radius_shrink > 0.15e-3)
                    {
                        sphere_radius_shrink -= 0.1 * sphere_radius;
                        next_direction = generateNextDirection(phi_to_target, theta_to_target, 0.03);
                        center         = next_direction * spheres_to_add[spheres_to_add.size() - 1].radius / 2.0 + last_sphere_center; 
                        if(isInVoxel(center, min_distance_from_border) && (!isSphereColliding(center, sphere_radius_shrink, center_to_test, last_radius)))
                        {
                            Sphere sphere_to_add(sphere_id, branch_id, center, sphere_radius_shrink);
                            if((sphere_to_add.center - spheres_to_add[spheres_to_add.size() - 1].center).norm() > 0.0004)
                            {
                                cout << "norm " << (sphere_to_add.center - spheres_to_add[spheres_to_add.size() - 1].center).norm() << endl;
                                cout << (center_to_test - center).norm() << endl;
                                cout << "CTT " << center_to_test << endl;
                                cout << "last " << last_sphere_center << endl;
                                cout << "s -2 " << spheres_to_add[spheres_to_add.size() - 2].center << endl;
                                cout << "s -1 " << spheres_to_add[spheres_to_add.size() - 1].center << endl;
                                cout << "s " << sphere_to_add.center << endl;
                            }

                            // cout << "p 3 " << sphere_to_add.center << endl;
                            spheres_to_add.push_back(sphere_to_add);
                            // cout << "p 3 " << spheres_to_add.back().center << endl;
                            // if(spheres_to_add.size() > 1)
                            //     cout << "p 3 " << (spheres_to_add[spheres_to_add.size() - 2].center - spheres_to_add[spheres_to_add.size() - 1].center).norm() << endl;

                            sphere_id++;
                            achieved = true;
                            last_sphere_center = spheres_to_add.back().center;
                            updateLUT(center, sphere_radius_shrink);
                            break;
                        }
                    }
                }
                
                
                // If we didn't manage to continue the branch but it's long enough
                double length_branch = (spheres_to_add[0].center - spheres_to_add[spheres_to_add.size() - 1].center).norm();
                if(!achieved && (length_branch > 15e-3))
                {
                    subbranch.set_spheres(spheres_to_add);
                    dendrite.add_subbranch(subbranch);
                    return {};
                }
                // Discard branch
                else if(!achieved && (length_branch < 15e-3))
                {
                    Sphere origin(0, 0, parent.origin, parent.radius);
                    spheres_to_add.insert(spheres_to_add.begin(), origin);
                    updateLUT(spheres_to_add, spheres_to_add.size() - 1);
                    return {};
                }
            } // isColliding
        }
        // Outside of the voxel. If long enough, keep the branch
        else
        {
            if(spheres_to_add.size() > 1)
            {
                double length_branch = (spheres_to_add[0].center - spheres_to_add[spheres_to_add.size() - 1].center).norm();
                if(length_branch > 30e-3)
                {
                    subbranch.set_spheres(spheres_to_add);
                    dendrite.add_subbranch(subbranch);
                }
                // Discard branch
                Sphere origin(0, 0, parent.origin, parent.radius);
                spheres_to_add.insert(spheres_to_add.begin(), origin);
                updateLUT(spheres_to_add, spheres_to_add.size() - 1);
                return {};
            }
            return {};
        } 
    }
    double length_branch = (spheres_to_add[0].center - spheres_to_add[spheres_to_add.size() - 1].center).norm();
    if(length_branch > 30e-3)
    {
        subbranch.set_spheres(spheres_to_add);
        dendrite.add_subbranch(subbranch);
    }

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
    
    // Return the next branching point
    return {spheres_to_add[spheres_to_add.size()-1].center, parent.direction, children_dir, branch_id, sphere_radius};
    
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

Vector3d NeuronDistribution::generatePointOnSphere() const
{
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::uniform_real_distribution<double> phi_dist(-1.0, 1.0);
    std::uniform_real_distribution<double> theta_dist(0.0, 2.0 * M_PI);

    double theta   = theta_dist(generator);
    double cos_phi = phi_dist(generator);
    double phi     = acos(cos_phi);
    double x       = sin(phi) * cos(theta);
    double y       = sin(phi) * sin(theta);
    double z       = cos(phi);

    return {x, y, z};   
}

Vector3d NeuronDistribution::generateNextDirection(double const& phi_to_target, double const& theta_to_target, double const& STD) const
{
    random_device rd{};
    mt19937 gen{rd()};

    // Generate tortuosity : theta ~ N(theta/(2*pi), STD)
    normal_distribution<double> theta_dist(theta_to_target / (2.0 * M_PI), STD);
    double theta = theta_dist(gen) * 2.0 * M_PI;
    if(theta > 2.0 * M_PI)
        theta = 2.0 * M_PI;
    else if (theta < -2.0 * M_PI)
        theta = -2.0 * M_PI;

    // Generate tortuosity : cos(phi) ~ N(cos(phi), STD*2*pi)
    double cos_phi_to_target = cos(phi_to_target);
    normal_distribution<double> cos_phi_dist(cos_phi_to_target, STD * 2.0 * M_PI);
    double cos_phi = cos_phi_dist(gen);
    if(cos_phi > 1)
        cos_phi = 1;
    else if (cos_phi < -1)
        cos_phi = -1;
    double phi = acos(cos_phi);
 
    double x     = sin(phi) * cos(theta);
    double y     = sin(phi) * sin(theta);
    double z     = cos(phi);
    
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

TEST_CASE("isSphereColliding")
{
    int num_obstacles = 1;
    double ICVF = 0;
    Vector3d min_limits_vx = {0.0, 0.0, 0.0};
    Vector3d max_limits_vx = {0.1, 0.1, 0.1};
    int sphere_overlap = 4;
    double step_length = 0.5e-3;
    NeuronDistribution neuronDist(num_obstacles, ICVF, min_limits_vx, max_limits_vx, sphere_overlap, step_length);

    // Add soma
    Vector3d soma_center = {0.05, 0.05, 0.05};
    double soma_radius   = 10e-3;
    neuronDist.updateLUT(soma_center, soma_radius);

    SUBCASE("non-colliding sphere")
    {
        Vector3d dir = {1.0, 0.0, 0.0};
        double dendrite_radius  = 0.5e-3;
        Vector3d center_to_test = soma_center + dir * (soma_radius + 3 * dendrite_radius);
        CHECK(!neuronDist.isSphereColliding(center_to_test, dendrite_radius));
    }
    SUBCASE("colliding sphere")
    {
        Vector3d dir = {1.0, 0.0, 0.0};
        double dendrite_radius  = 0.5e-3;
        Vector3d center_to_test = soma_center + dir * soma_radius;
        CHECK(neuronDist.isSphereColliding(center_to_test, dendrite_radius));
    }
    
}

bool NeuronDistribution::isSphereColliding(Vector3d const& sphere_center, double const& sphere_radius) 
{
    uint16_t min_x, min_y, min_z, max_x, max_y, max_z;
    min_x = floor((sphere_center[0] - sphere_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    max_x = floor((sphere_center[0] + sphere_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    min_y = floor((sphere_center[1] - sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    max_y = floor((sphere_center[1] + sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    min_z = floor((sphere_center[2] - sphere_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
    max_z = floor((sphere_center[2] + sphere_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
    // cout << min_x << " " << max_x << " " << min_y << " " << max_y << " " << min_z << " " << max_z << endl;

    for(size_t x=min_x; x < max_x; ++x)
        for(size_t y=min_y; y < max_y; ++y)
            for(size_t z=min_z; z < max_z; ++z)
                if(LUT[x][y][z])
                    return true;
    return false; 
}

TEST_CASE("isSphereColliding_ignore")
{
    int num_obstacles = 1;
    double ICVF = 0;
    Vector3d min_limits_vx = {0.0, 0.0, 0.0};
    Vector3d max_limits_vx = {0.1, 0.1, 0.1};
    int sphere_overlap = 4;
    double step_length = 0.5e-3;
    NeuronDistribution neuronDist(num_obstacles, ICVF, min_limits_vx, max_limits_vx, sphere_overlap, step_length);

    // Add soma
    Vector3d soma_center = {0.05, 0.05, 0.05};
    double soma_radius   = 10e-3;
    neuronDist.updateLUT(soma_center, soma_radius);

    SUBCASE("non-colliding sphere (non-touching)")
    {
        Vector3d dir = {1.0, 0.0, 0.0};
        double dendrite_radius  = 0.5e-3;
        Vector3d center_to_test = soma_center + dir * (soma_radius + 3 * dendrite_radius);
        CHECK(!neuronDist.isSphereColliding(center_to_test, dendrite_radius, soma_center, soma_radius));
    }
    SUBCASE("non-colliding sphere (touching)")
    {
        Vector3d dir = {1.0, 0.0, 0.0};
        double dendrite_radius  = 0.5e-3;
        Vector3d center_to_test = soma_center + dir * soma_radius;
        CHECK(!neuronDist.isSphereColliding(center_to_test, dendrite_radius, soma_center, soma_radius));
    }
    SUBCASE("dendrite")
    {
        Vector3d dendrite_center;
        Vector3d dir = {1.0, 0.0, 0.0};
        double dendrite_radius  = 0.5e-3;
        Vector3d last_center = soma_center + dir * soma_radius;
        for(size_t i=1; i < 50; ++i)
        {
            dendrite_center = soma_center + dir * (soma_radius + double(i) * dendrite_radius / 2.0);
            CHECK(!neuronDist.isSphereColliding(dendrite_center, dendrite_radius, last_center, dendrite_radius));
            neuronDist.updateLUT(dendrite_center, dendrite_radius);
            last_center = dendrite_center;
        }
    }   
}

bool NeuronDistribution::isSphereColliding(Vector3d const& sphere_center, double const& sphere_radius, Vector3d const& sphere_to_ignore_center, double const& sphere_to_ignore_radius) 
{
    uint16_t min_x = floor((sphere_center[0] - sphere_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    uint16_t max_x = floor((sphere_center[0] + sphere_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    uint16_t min_y = floor((sphere_center[1] - sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    uint16_t max_y = floor((sphere_center[1] + sphere_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    uint16_t min_z = floor((sphere_center[2] - sphere_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
    uint16_t max_z = floor((sphere_center[2] + sphere_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));


    uint16_t min_ign_x = floor((sphere_to_ignore_center[0] - sphere_to_ignore_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    uint16_t max_ign_x = floor((sphere_to_ignore_center[0] + sphere_to_ignore_radius) / ((max_limits_vx[0] - min_limits_vx[0]) / SIZE_LUT));
    uint16_t min_ign_y = floor((sphere_to_ignore_center[1] - sphere_to_ignore_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    uint16_t max_ign_y = floor((sphere_to_ignore_center[1] + sphere_to_ignore_radius) / ((max_limits_vx[1] - min_limits_vx[1]) / SIZE_LUT));
    uint16_t min_ign_z = floor((sphere_to_ignore_center[2] - sphere_to_ignore_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
    uint16_t max_ign_z = floor((sphere_to_ignore_center[2] + sphere_to_ignore_radius) / ((max_limits_vx[2] - min_limits_vx[2]) / SIZE_LUT));
    
    // cout << min_x << " " << max_x << " " << min_y << " " << max_y << " " << min_z << " " << max_z << endl;
    // cout << min_ign_x << " " << max_ign_x << " " << min_ign_y << " " << max_ign_y << " " << min_ign_z << " " << max_ign_z << endl;
    for(size_t x=min_x; x < max_x; ++x)
        for(size_t y=min_y; y < max_y; ++y)
            for(size_t z=min_z; z < max_z; ++z)
                if(LUT[x][y][z])
                {
                    if(!((min_ign_x <= x) && (max_ign_x >= x) &&
                         (min_ign_y <= y) && (max_ign_y >= y) &&
                         (min_ign_z <= z) && (max_ign_z >= z)))
                         {
                            return true;      
                         }
                }
                    
    return false; 
}

bool NeuronDistribution::isSphereCollidingSphere(Vector3d const& pos1, Vector3d const& pos2, double const& radius1, double const& radius2, double const& minDistance) const 
{
    Vector3d m = pos1 - pos2;
    double distance_to_sphere = m.norm() - radius1 - radius2;

    return distance_to_sphere < minDistance;
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
        out << "Neuron " + to_string(i) << endl;
        out << "Soma " + to_string(i) << endl;
        // Print for soma : x y z r bool_active
        out << neurons[i].soma.center[0] << " " 
        << neurons[i].soma.center[1] << " "
        << neurons[i].soma.center[2] << " "
        << neurons[i].soma.radius    << endl; 

        for (size_t j = 0; j < neurons[i].dendrites.size(); j++)
        {
            out << "Dendrite " + to_string(j) << endl;
            for (size_t k = 0; k < neurons[i].dendrites[j].subbranches.size(); k++)
            {
                out << "Segment " + to_string(k) << endl;
                for (size_t l = 0; l < neurons[i].dendrites[j].subbranches[k].spheres.size(); l++)
                {
                    // Print for each dendrite, each sphere
                    out << neurons[i].dendrites[j].subbranches[k].spheres[l].center[0] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[1] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].center[2] << " "
                    << neurons[i].dendrites[j].subbranches[k].spheres[l].radius << endl; 
                        
                }
                 
                if(neurons[i].dendrites[j].subbranches[k].proximal_branching.size() == 2)
                    out << "proximal " + to_string(neurons[i].dendrites[j].subbranches[k].proximal_branching[0]) + " " + to_string(neurons[i].dendrites[j].subbranches[k].proximal_branching[1]);
                else
                    out << "proximal " + to_string(neurons[i].dendrites[j].subbranches[k].proximal_branching[0]);

                if(neurons[i].dendrites[j].subbranches[k].distal_branching.size() == 2)
                    out << " distal " + to_string(neurons[i].dendrites[j].subbranches[k].distal_branching[0]) + " " + to_string(neurons[i].dendrites[j].subbranches[k].distal_branching[1]) << endl;
                else
                    out << " distal -1" << endl;
            }
        }
        out << "end" << endl;
    }
}

/* Function that orders the indices of the tree from the origin to the end
   for all end nodes. */
vector<int> order_tree_idx(vector<vector<int>> const& paths) 
{
    vector<int> order;
    size_t idx_longest_branch = 0;
    size_t max_branch_length  = 0;
    for(size_t i=0; i < paths.size(); ++i)
    {
        if(paths[i].size() > max_branch_length)
        {
            max_branch_length  = paths[i].size();
            idx_longest_branch = i;
        }
    }
    for(size_t i=1; i < paths[idx_longest_branch].size(); ++i)
    {
        order.push_back(paths[idx_longest_branch][i]);
    }
    order.push_back(paths[idx_longest_branch + 1][paths[idx_longest_branch + 1].size() - 1]);
    order.push_back(paths[idx_longest_branch + 2][paths[idx_longest_branch + 2].size() - 2]);
    order.push_back(paths[idx_longest_branch + 2][paths[idx_longest_branch + 2].size() - 1]);
    order.push_back(paths[idx_longest_branch + 3][paths[idx_longest_branch + 3].size() - 1]);

    return order;
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

    out << 1 << endl; //scale
    out << 0 << endl; //volume_inc_perc
    out << 0 << endl; //dyn_perc
    out << icvf << endl;
    out << min_limits_vx[0] << endl; //min_limits [mm]
    out << max_limits_vx[0] << endl; //max_limits [mm]

    for (size_t i = 0; i < neurons.size(); i++)
    {
        int line  = 1;
    
        // SampleID TypeID x y z r ParentID
        out << line                 << " " 
        << 1                         << " "
        << neurons[i].soma.center[0] << " "
        << neurons[i].soma.center[1] << " "
        << neurons[i].soma.center[2] << " "
        << neurons[i].soma.radius    << " "
        << -1                        << endl; 
        line++;

    
        // TODO [ines] : remove this ugly hard code...
        vector<int> line_node = {3, 4, 4, 3, 7, 7};
        for(size_t dendrite_id=0; dendrite_id < neurons[i].dendrites.size(); ++dendrite_id)
        {
            auto first_subbranch = (neurons[i]).dendrites[dendrite_id].subbranches[0];
            auto sphere_begin    = first_subbranch.spheres[0];
            auto sphere_end      = first_subbranch.spheres[first_subbranch.spheres.size() - 1];
            
            // Soma to first node, at the soma surface
            out << line                 << " " 
            << 3                        << " "
            << sphere_begin.center[0]   << " "
            << sphere_begin.center[1]   << " "
            << sphere_begin.center[2]   << " "
            << sphere_begin.radius      << " "
            << 1                        << endl; 
            line++;

            // First subbranch
            out << line                 << " " 
            << 3                        << " "
            << sphere_end.center[0]     << " "
            << sphere_end.center[1]     << " "
            << sphere_end.center[2]     << " "
            << sphere_end.radius        << " "
            << line - 1                 << endl; 
            line++;

            vector<vector<int>> paths = find_tree_paths(new Neuron(neurons[i]), dendrite_id);
            vector<int> order = order_tree_idx(paths);

            for(size_t subbranch_id=0; subbranch_id < order.size(); ++subbranch_id)
            {
                auto spheres     = (neurons[i]).dendrites[dendrite_id].subbranches[order[subbranch_id]].spheres;
                auto sphere_node = spheres[spheres.size() - 1];
                
                // SampleID TypeID x y z r ParentID
                out << line                 << " " 
                << 3                        << " "
                << sphere_node.center[0]    << " "
                << sphere_node.center[1]    << " "
                << sphere_node.center[2]    << " "
                << sphere_node.radius       << " "
                << line_node[subbranch_id] + dendrite_id * 8 << endl; 
                line++;
            }
        }
    }
}


/* Find all the paths from origin to end node */
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

/* Find the path from origin to end node */
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

bool intersection_sphere_vector(double const& min_distance, Vector3d const& sphere_center, double const& sphere_radius, Vector3d const& step_dir, Vector3d const& traj_origin)
{
    Vector3d m = traj_origin - sphere_center;
    double rad = sphere_radius;

    int a = 1;
    double b = m.dot(step_dir);
    double c = m.dot(m) - rad * rad;

    double discr = b * b - a * c;

    if (discr < 0.0)
        return false;

    auto intercept1 = (-b + sqrt(discr)) / a;
    auto intercept2 = (-b - sqrt(discr)) / a;

    if((intercept1 < min_distance && intercept1 > 0) || 
       (intercept2 < min_distance && intercept2 > 0))
        return true;
    else
        return false;
}

bool NeuronDistribution::avoid_somas(double const& min_distance, vector<Vector3d> soma_centers, vector<double> const& soma_radii, Vector3d const& soma_to_ignore, Vector3d const &step_dir, Vector3d const &traj_origin) const
{
    for(size_t soma_id=0; soma_id < soma_centers.size(); ++soma_id)
    {
        if(soma_radii[soma_id] == 10e-3)
            if(((soma_centers[soma_id]-soma_to_ignore).norm() > 2.0 * soma_radii[soma_id]) &&
                intersection_sphere_vector(min_distance, soma_centers[soma_id], soma_radii[soma_id], step_dir, traj_origin))
                return false;
    }
    return true;
}

