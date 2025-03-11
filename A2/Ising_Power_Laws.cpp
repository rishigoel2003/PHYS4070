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
    int L = 128;
    double B=0;
    int N = pow(L,2);

    
    vector<vector<int>> lattice(L, vector<int>(L));

    initializeLattice(lattice);

    int numTemperatures = 30;
    double minTemp = 2.1; // 
    double maxTemp = 2.26;


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
        
        
        int sweeps = 10000; // Number of iterations
        



        // Run simulation
        for (int i = 0; i < sweeps*pow(L,2); i++) {
            
        


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
            

            if (i >= (3*sweeps/4)*N && i % 10*N == 0) {
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





