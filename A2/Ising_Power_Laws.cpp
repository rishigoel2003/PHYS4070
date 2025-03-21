#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>

std::ofstream output_file("IsingPower.txt");

using namespace std;

#include <random>

#include "initialise.hpp"
#include "functions.hpp"



// Main function
int main() {
    // large lattice size required for this part
    int L = 256;
    double B=0;
    int N = pow(L,2);

    
    vector<vector<int>> lattice(L, vector<int>(L));

    initializeLattice(lattice);

    //lots of temperatures required to reproduce power laws
    int numTemperatures = 100;
    double minTemp = 2; // 
    double maxTemp = 2.27;


    double energy = ComputeInitialEnergy(lattice);
    double magnetisation = computeInitialMagnetization(lattice);
    cout << "Initial energy: " << energy << " Initial Magnetisation: " << magnetisation << endl; 

    // printLattice(lattice);

    random_device rd;
    mt19937 gen(rd());
    
    uniform_int_distribution<> distrib(0, L-1);
    uniform_real_distribution<double> dist(0.0, 1.0);


    for (int tempIndex = 0; tempIndex < numTemperatures; tempIndex++) {
        // Calculate current temperature
        double Temperature = minTemp + tempIndex * (maxTemp - minTemp) / (numTemperatures - 1);
        
        cout << "Running simulation at T = " << Temperature << endl;
        
        
        int sweeps = 2000; // Number of iterations
        



        // Run simulation
        for (int i = 0; i < sweeps*N; i++) {
            
        


            // Choose a random position in the lattice
            
            int randomRow = distrib(gen);
            int randomCol = distrib(gen);
            
            // Calculate energy change if this spin is flipped
            double deltaE = calculateDeltaE(lattice, randomRow, randomCol,B);
            
            // Decide whether to flip the spin based on the Metropolis algorithm
            if (deltaE <= 0 || (exp(-deltaE/Temperature) >= dist(gen))) {
                // Flip the spin
                lattice[randomRow][randomCol] *= -1;
                // Update energy
                energy += deltaE;
                //update magnetisation
                magnetisation += 2.0*lattice[randomRow][randomCol];
            }
            

            if (i >= (5*sweeps/10)*N && i % (5*N) == 0) {
                // Take measurements here
                output_file << i/N << " " << Temperature << " " << energy/N << " " << magnetisation/N << "\n";
            }

    }


    // cout << "algorithm finished\n";

    // printLattice(lattice);

    // cout << "Final energy: " << energy << " Final Magnetisation: " << magnetisation << endl; 
    
    }
    

    
    output_file.close();

    return 0;
}





