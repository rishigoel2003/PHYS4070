#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>

// name of output file
std::ofstream output_file("trajectories.txt");


using namespace std;


//including the helper functions and class
#include "particle_class.hpp"
#include "functions.hpp"






// Main function
int main() {
    //setting constants found analytically
    const double moon_initial_speed = sqrt(1/19.0);
    const double projectile_initial_speed = sqrt(142/75.0);
    const double moon_orbit_radius = 19;
    const double rads_below_x = 1.0798828009;

    Particle planet = {{0, 0}, {0, 0}, 1.0}; // Planet at origin with mass 1

    // Moon at starting position with direction for exact meeting with projectile, mass 0.1
    Particle moon = {{cos(rads_below_x)*moon_orbit_radius,-sin(rads_below_x) * moon_orbit_radius}, {sin(rads_below_x) * moon_initial_speed,cos(rads_below_x)*moon_initial_speed}, 0.1};

    // Projectile starts on planet right side and moves radially outward along x axis with negligible mass
    Particle projectile = {{1, 0}, {projectile_initial_speed,0}, 0};
    
    
    double dt = 0.0001;  // Time step
    int steps = 1200000; // Number of iterations
    double moon_proj_distance = sqrt(pow((moon.r[0]- projectile.r[0]),2) + pow((moon.r[1]- projectile.r[1]),2)); //calculate euclidean distance

    // Run simulation
    for (int i = 0; i < steps; i++) {
        //stepping forward both moon and planet as they are independent
        rk4_step(moon, planet, dt);
        rk4_step(projectile, planet, dt);

        //update distance
        moon_proj_distance = sqrt(pow((moon.r[0]- projectile.r[0]),2) + pow((moon.r[1]- projectile.r[1]),2));

        //if we have done enough steps print outputs to file
        if (i%10000 == 1){
        // Output the positions to the file
        output_file << i * dt << " " << moon.r[0] << " " << moon.r[1] << " "
                    << projectile.r[0] << " " << projectile.r[1] << " " << moon_proj_distance << "\n";
                    }

        //checking if the projectile is at the surface of the moon (or gets close enough due to numerical rounding)
        if (projectile.r[0] >= 18.7499) {
            cout << "Time when projectile makes it to the moon: " << i*dt << "\n";
            cout << "position of the moon: (" << moon.r[0] << "," << moon.r[1] << ")" << "\n"; 
            break;
        }
    }
    
    //close output file
    output_file.close();

    // end main
    return 0;
}



