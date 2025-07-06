#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>

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
    const double tMax = 50.0;           // Maximum simulation time
    const int numSteps = tMax / dt;     // Number of time steps
    const int saveInterval = 20;       // Save data every saveInterval steps
    
    vector<double> x = initialise_grid(N, L, dx);
    
    // Initialize wavefunction as empty vector
    ComplexVector psi(N, 0.0);
    
    double g = -1;              // Interaction strength (attractive)
    double u = 0.1;           // u values


    vector<double> theta_vals = {2*pi,pi,2.0*pi/3.0,pi/2.0,0}; //theta values for relative phase

    ofstream outFile("Part1/nlse_evolution_solitons.dat");
    ofstream outFile2("Part1/nlse_evolution_solitons_peaks.dat");


    //loop over values of u
    for (int i = 0; i < theta_vals.size(); ++i) {
        double theta = theta_vals[i];


                
        //initialise wavefunction with a plane wave
        solitons(psi,x,u,theta,N,L,dx);

            // Save initial state
        for (int i = 0; i < N; ++i) { // pre incrementing is meant to be good practise
            outFile << theta << " " <<  x[i] << " " << 0 << " " << norm(psi[i]) << endl;
        }

        outFile << endl;
        
            
        
        // Time evolution
        cout << "theta=" << theta << endl;
        for (int step = 1; step <= numSteps; ++step) {
            // Apply RK4 step
            psi = rungeKutta4(psi, dt, dx, g);
        
            // Output data at certain intervals
            if (step % saveInterval == 0) {
                double currentTime = step * dt;
        
                // Inside your loop
                std::vector<std::pair<double, int>> localPeaks; // (density, index)

                for (int i = 1; i < N - 1; ++i) {
                    double left = norm(psi[i - 1]);
                    double center = norm(psi[i]);
                    double right = norm(psi[i + 1]);

                    if (center > left && center > right) {
                        localPeaks.emplace_back(center, i);
                    }

                    // Save wavefunction (can stay as is)
                    outFile << theta << " " << x[i] << " " << currentTime << " " << center << endl;
                }
                outFile << endl;

                // Sort peaks by descending probability density
                sort(localPeaks.begin(), localPeaks.end(), std::greater<>());

                if (localPeaks.size() >= 2) {
                    int maxIndex1 = localPeaks[0].second;
                    int maxIndex2 = localPeaks[1].second;
                    double peakDistance = std::abs(x[maxIndex1] - x[maxIndex2]);

                    outFile2 << theta << " " << currentTime << " "
                            << min(x[maxIndex1],x[maxIndex2]) << " " << max(x[maxIndex2],x[maxIndex1]) << " " << peakDistance << std::endl;
                }
            }
        }
        


    }

outFile.close();
cout << "Simulation completed." << endl;

return 0;
        
    

    
}