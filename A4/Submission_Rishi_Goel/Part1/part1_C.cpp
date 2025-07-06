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
    const double tMax = 30.0;           // Maximum simulation time
    const int numSteps = tMax / dt;     // Number of time steps
    const int saveInterval = 10;       // Save data every saveInterval steps
    
    vector<double> x = initialise_grid(N, L, dx);
    
    // Initialize wavefunction as empty vector
    ComplexVector psi(N, 0.0);
    
    double g = -1;              // Interaction strength (attractive)

    //values for g
    vector<double> u_vals = {-1,-0.5, 0,0.5,1}; 
    // vector<double> u_vals = {-1,0,1}; // 


    ofstream outFile("Part1/nlse_evolution_wave_packet.dat");
    if (!outFile) {
        cerr << "Failed to open output file." << endl;
        return 1;
    }

    ofstream outFile2("Part1/nlse_evolution_wave_packet_peaks.dat");


    //loop over values of u
    for (int i = 0; i < u_vals.size(); ++i) {
        double u = u_vals[i];


                
        //initialise wavefunction with a plane wave
        ComplexVector psi(N, 0.0);  // fresh psi
        wave_packet(psi, x, u, N, dx);


            // Save initial state
        for (int i = 0; i < N; ++i) { // pre incrementing is meant to be good practise
            outFile << u << " " <<  x[i] << " " << 0 << " " << norm(psi[i]) << endl;
        }

        outFile << endl;
        
            
        
        // Time evolution
        cout << "u=" << u << endl;
        for (int step = 1; step <= numSteps; ++step) {
            // Apply RK4 step
            psi = rungeKutta4(psi, dt, dx, g);
        
            // Output data at certain intervals
            if (step % saveInterval == 0) {
                double currentTime = step * dt;
        
                double maxProbDensity = 0.0;
                int maxIndex = 0;
        
                // int maxV = max_element(V.begin(), V.end())−V.begin()

                for (int i = 0; i < N; ++i) {
                    //finding the maximum probability density
                    double probDensity = norm(psi[i]); // |psi|^2
                    if (probDensity > maxProbDensity) {
                        maxProbDensity = probDensity;
                        maxIndex = i;
                    }
        
                    // Save wavefunction
                    outFile << u << " " << x[i] << " " << currentTime << " " << probDensity << endl;
                }
                outFile << endl;
        
                // Save peak position
                outFile2 << u << " " << currentTime << " " << x[maxIndex] << endl;
            }
        }
        


    }

outFile.close();
cout << "Simulation completed." << endl;

return 0;
        
    

    
}