#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>

std::ofstream output_file("trajectories_B2.txt");


using namespace std;

#include "particle_class.hpp"
#include "functions.hpp"


// Main function
int main() {
    //same initial conditions as B1
    const double moon_initial_speed = sqrt(1/19.0);
    const double projectile_initial_speed = sqrt(142/75.0);
    const double moon_orbit_radius = 19;
    const double rads_below_x = 1.0798828009;

    Particle planet = {{0, 0}, {0, 0}, 1.0}; // Planet at origin with mass 1

    // initial conditions for moon
    Particle moon = {{cos(rads_below_x)*moon_orbit_radius,-sin(rads_below_x) * moon_orbit_radius}, {sin(rads_below_x) * moon_initial_speed,cos(rads_below_x)*moon_initial_speed}, 0.1}; 

    // Projectile starts on surface and moves radially outward, mass 0 as it is negligible (within floating point precision error)
    Particle projectile = {{1, 0}, {projectile_initial_speed,0}, 0}; 

    double dt = 0.0001;  // Time step
    int steps = 1200000; // Number of iterations

    //calculate distance
    double moon_proj_distance = sqrt(pow((moon.r[0]- projectile.r[0]),2) + pow((moon.r[1]- projectile.r[1]),2));


    // Run simulation
    for (int i = 0; i < steps; i++) {
        //update both projectile and moon simultaneously
        rk4_step_3body(projectile, moon,planet, dt);
        
        //update distance
        moon_proj_distance = sqrt(pow((moon.r[0]- projectile.r[0]),2) + pow((moon.r[1]- projectile.r[1]),2));
        const double epsilon = 5e-5; // Small tolerance

        //check if the moon is within some small tolerance of 0.5 radius
        if (fabs(moon_proj_distance - 0.5) < epsilon) {
            cout << "Time when projectile is 0.5 units away from moon:" << i * dt << "\n";
        }



        if (i%10000 == 1){
        // Output the positions to the file
        output_file << i * dt << " " << moon.r[0] << " " << moon.r[1] << " "
                    << projectile.r[0] << " " << projectile.r[1]<< " " << moon_proj_distance << "\n";
                    }

    }
    
    output_file.close();


    return 0;
}



