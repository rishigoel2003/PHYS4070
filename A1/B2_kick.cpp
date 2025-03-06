#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <math.h>

std::ofstream output_file("trajectories_B2_kick.txt");


using namespace std;

#include "particle_class.hpp"
#include "functions.hpp"





// Main function
int main() {

    // same initial conditions and B1 and B2
    const double moon_initial_speed = sqrt(1/19.0);
    const double projectile_initial_speed = sqrt(142/75.0);
    const double moon_orbit_radius = 19;
    const double rads_below_x = 1.0798828009;
    
    Particle planet = {{0, 0}, {0, 0}, 1.0}; // Planet at origin with mass 1
    Particle moon = {{cos(rads_below_x)*moon_orbit_radius,-sin(rads_below_x) * moon_orbit_radius}, {sin(rads_below_x) * moon_initial_speed,cos(rads_below_x)*moon_initial_speed}, 0.1}; 
    Particle projectile = {{1, 0}, {projectile_initial_speed,0}, 0}; 


    double dt = 0.0001;  // Time step
    int steps = 4000000; // Number of iterations

    double moon_proj_distance = 0;



    double epsilon = 1e-4; // Small tolerance

    // Run simulation
    for (int i = 0; i < steps; i++) {
        rk4_step_3body(projectile, moon,planet, dt); //update position of both projectile and moon
        
        //update distance
        moon_proj_distance = sqrt(pow((moon.r[0]- projectile.r[0]),2) + pow((moon.r[1]- projectile.r[1]),2));
        

        // if within the tolerance of 0.5 radius, apply kick
        if (fabs(moon_proj_distance - 0.5) < epsilon) {
            
            //calculate angular distance
            double theta = fabs(atan2(moon.r[1] - projectile.r[1], moon.r[0] - projectile.r[0]));
            
            //define the kick velocity based on analytical calculation
            double kick_x = moon.v[0] + 1/sqrt(5)*sin(theta) - projectile.v[0];
            double kick_y = moon.v[1] + 1/sqrt(5)*cos(theta) - projectile.v[1];

            //apply kick
            projectile.v[0] += kick_x;
            projectile.v[1] += kick_y;
            
            //set epsilon negative so this condition doesnt trigger again
            epsilon = -1;
            cout << "Projectile is 0.5 away from moon, kick activated for orbital transfer. \n";
        }

        if (i%5000 == 1){
        // Output the positions to the file
        output_file << i * dt << " " << moon.r[0] << " " << moon.r[1] << " "
                    << projectile.r[0] << " " << projectile.r[1]<< " " << moon_proj_distance << "\n";
                    }

    }
    
    output_file.close();

    return 0;
}



