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
    const double moon_initial_speed = sqrt(1/19.0);
    const double projectile_initial_speed = sqrt(142/75.0);
    const double moon_orbit_radius = 19;
    const double rads_below_x = 1.0798828009;

    Particle planet = {{0, 0}, {0, 0}, 1.0}; // Planet at origin with mass 1

    Particle moon = {{cos(rads_below_x)*moon_orbit_radius,-sin(rads_below_x) * moon_orbit_radius}, {sin(rads_below_x) * moon_initial_speed,cos(rads_below_x)*moon_initial_speed}, 0.1}; // Moon at (19,0) with velocity in y-direction


    Particle projectile = {{1, 0}, {projectile_initial_speed,0}, 0}; // Projectile starts on surface and moves radially outward

    double dt = 0.0001;  // Time step
    int steps = 1200000; // Number of iterations

    double moon_proj_distance = sqrt(pow((moon.r[0]- projectile.r[0]),2) + pow((moon.r[1]- projectile.r[1]),2));
;

    // Run simulation
    for (int i = 0; i < steps; i++) {
        rk4_step_3body(projectile, planet, moon, dt);
        rk4_step_3body(moon, planet, projectile, dt);
        
        moon_proj_distance = sqrt(pow((moon.r[0]- projectile.r[0]),2) + pow((moon.r[1]- projectile.r[1]),2));
        const double epsilon = 1e-4; // Small tolerance

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



