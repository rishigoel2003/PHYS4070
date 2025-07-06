#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
using Complex = complex<double>;
using ComplexVector = vector<Complex>;

#include "functions.hpp"




int main() {
    // Parameters
    const double L = 20.0;              // Domain size
    const int N = 128;                  // Number of grid points
    const double dx = L / (N - 1);      // Spatial step size
    const double dt = 0.01;             // Time step
    const double tMax = 20.0;           // Maximum simulation time
    const int numSteps = tMax / dt;     // Number of time steps
    const int saveInterval = 10;       // Save data every saveInterval steps
    
    vector<double> x = initialise_grid(N, L, dx);
    
    // Initialize wavefunction as empty vector
    ComplexVector psi(N, 0.0);
    
    //values for g
    vector<double> g_vals = {-10,-3,-2,-1,-0.5, 0}; // Interaction strength (attractive)


    ofstream outFile("Part1/nlse_evolution_plane_wave.dat");
    if (!outFile) {
        cerr << "Failed to open output file." << endl;
        return 1;
    }


    //loop over values of g
    for (int i = 0; i < g_vals.size(); ++i) {
        double g = g_vals[i];


                
        //initialise wavefunction with a plane wave
        plane_wave(psi,x,N,L,dx);

            // Save initial state
        for (int i = 0; i < N; ++i) { // pre incrementing is meant to be good practise
            outFile << g << " " <<  x[i] << " " << 0 << " " << norm(psi[i]) << endl;
        }

        outFile << endl;
        
            
        
        // Time evolution
        cout << "g=" << g << endl;
        for (int step = 1; step <= numSteps; ++step) {
            // Apply RK4 step
            psi = rungeKutta4(psi, dt, dx, g);
            
            // Output data at certain intervals
            if (step % saveInterval == 0) {
                double currentTime = step * dt;
                // double currentNorm = calculateNorm(psi, dx);
                
                    // Save wavefunction
                for (int i = 0; i < N; ++i) {
                    outFile << g << " " <<  x[i] << " " << currentTime << " " << norm(psi[i]) << endl;
                }
                outFile << endl;
                
            }
        }


    }

outFile.close();
cout << "Simulation completed." << endl;

return 0;
        
    

    
}