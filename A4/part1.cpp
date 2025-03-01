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
    const double tMax = 10.0;           // Maximum simulation time
    const int numSteps = tMax / dt;     // Number of time steps
    const double g = -1.0;              // Interaction strength (attractive)
    const int saveInterval = 100;       // Save data every saveInterval steps
    
    vector<double> x = initialise_grid(N, L, dx);
    
    // Initialize wavefunction with a Gaussian pulse
    ComplexVector psi(N, 0.0);
    
    double A = 1;                     // Amplitude
    double sigma = 1;                 // Width
    double x0 = 0;                    // Center
    double k0 = 0;                    // Initial momentum
    
    gaussian(psi, x, A, sigma, x0, k0, N, dx);

    
    // Output initial state
    ofstream outFile("nlse_evolution.dat");
    if (!outFile) {
        cerr << "Failed to open output file." << endl;
        return 1;
    }
    
    // Save initial state
    for (int i = 0; i < N; ++i) { // pre incrementing is meant to be good practise
        outFile << x[i] << " " << 0.0 << " " << abs(psi[i]) << " " 
                << real(psi[i]) << " " << imag(psi[i]) << endl;
    }

    
    outFile << endl;
    
    // Time evolution
    cout << "Starting simulation..." << endl;
    for (int step = 1; step <= numSteps; ++step) {
        // Apply RK4 step
        psi = rungeKutta4(psi, dt, dx, g);
        
        // Output data at certain intervals
        if (step % saveInterval == 0) {
            double currentTime = step * dt;
            double currentNorm = calculateNorm(psi, dx);
            
            cout << "Time: " << currentTime << ", Norm: " << currentNorm << endl;
            
            // Save wavefunction
            for (int i = 0; i < N; ++i) {
                outFile << x[i] << " " << currentTime << " " << abs(psi[i]) << " " 
                        << real(psi[i]) << " " << imag(psi[i]) << endl;
            }
            outFile << endl;
        }
    }
    
    outFile.close();
    cout << "Simulation completed. Results saved to 'nlse_evolution.dat'." << endl;
    
    return 0;
}