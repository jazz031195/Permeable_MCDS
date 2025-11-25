// main_traj.cpp
#include <string>

int run_trajectory(const std::string& trajectoryFile) {
    // Initialize simulation
    Simulation sim;
    sim.loadConfiguration("config.conf");
    sim.initialize();

    // Run simulation
    sim.run();

    // Write trajectory to file
    sim.writeTrajectory(trajectoryFile);

    return 0;
}