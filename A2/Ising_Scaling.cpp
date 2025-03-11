#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>

std::ofstream output_file("IsingScaling.txt");

using namespace std;

#include <random>

#include "initialise.hpp"
#include "functions.hpp"



// Main function
int main() {

    random_device randDevice;
    mt19937 gen(randDevice());

    double B=0; //magnetic field

    vector<int> sizes = {8,16,32,64};

    int number_sizes = size(sizes);

    for (int l=0;l<number_sizes;l++){
    
        int L = sizes[l];
        int N = pow(L,2);

        
        vector<vector<int>> lattice(L, vector<int>(L));

        initializeLattice(lattice);

        int numTemperatures = 30;
        double minTemp = 0.3; // Avoiding T=0 for stability
        double maxTemp = 5.0;


        double energy = ComputeInitialEnergy(lattice);
        double magnetisation = computeInitialMagnetization(lattice);
        cout << "Initial energy: " << energy << " Initial Magnetisation: " << magnetisation << endl; 

        // printLattice(lattice);

        
        


        for (int tempIndex = 0; tempIndex < numTemperatures; tempIndex++) {
            // Calculate current temperature
            double Temperature = minTemp + tempIndex * (maxTemp - minTemp) / (numTemperatures - 1);
            
            cout << "Running simulation at T = " << Temperature << ", L = " << L << endl;
            
            
            int sweeps = 2000; // Number of iterations
            

            // Run simulation
            for (int i = 0; i < sweeps*N; i++) {
                
            


                // Choose a random position in the lattice
                
                uniform_int_distribution<> distrib(0, L-1);
                int randomRow = distrib(gen);
                int randomCol = distrib(gen);
                
                // Calculate energy change if this spin is flipped
                double deltaE = calculateDeltaE(lattice, randomRow, randomCol, B);
                
                
                uniform_real_distribution<double> dist(0.0, 1.0);
                // Decide whether to flip the spin based on the Metropolis algorithm
                if (deltaE <= 0 || (exp(-deltaE/Temperature) >= dist(gen))) {
                    // Flip the spin
                    lattice[randomRow][randomCol] *= -1;
                    // Update energy
                    energy += deltaE;
                    //update magnetisation
                    magnetisation += 2.0*lattice[randomRow][randomCol];
                }
                

                if (i >= (sweeps/2)*N && i % (L*L) == 0)//condition found from part 1 (1000 sweeps for burn in and 1000 for average)
                {
                    // Take measurements here
                    output_file << i/(N) << " "<< L << " " << Temperature << " " << energy/N << " " << magnetisation/N << "\n";
                }

        }


        // cout << "algorithm finished\n";

        printLattice(lattice);

        cout << "Final energy: " << energy << " Final Magnetisation: " << magnetisation << endl; 
        
        }
        
    }
    
    output_file.close();

    return 0;
}


