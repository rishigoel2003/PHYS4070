#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>
#include <omp.h>  // Include OpenMP header

using namespace std;

#include <random>

#include "initialise.hpp"
#include "functions.hpp"

// Main function
int main() {
    double t1 = omp_get_wtime();

    int L = 64;
    
    vector<vector<int>> lattice(L, vector<int>(L));

    initializeLattice(lattice);
    double energy = ComputeInitialEnergy(lattice);
    double magnetisation = computeInitialMagnetization(lattice);
    
    int sweeps = 10000; // Number of iterations
    double Temperature = 5;


    

    // Run simulation with OpenMP
    #pragma omp parallel default(none) shared(energy, magnetisation, Temperature, sweeps, lattice, L, std::cout)
    {

        int thread_id = omp_get_thread_num();

        // Generate a "thread local" random number generator
        thread_local std::mt19937 generator_tl(std::random_device{}() +
        omp_get_thread_num());

        uniform_real_distribution<double> prob_dist = uniform_real_distribution<double>(0.0, 1.0);
        
        // Each thread processes a portion of the total iterations
        for (int i = 0; i < sweeps; i++) {
            // Choose a random position in the lattice
            double dE = 0;
            double dM = 0;
            
            
            #pragma omp for 
            for (int row = 0; row<L; row+=2){
                for (int column = 0; column<L; column++){
                // Calculate energy change if this spin is flipped
                double deltaE = calculateDeltaE(lattice, row, column);
                
                // Decide whether to flip the spin based on the Metropolis algorithm
                if (deltaE <= 0 || (exp(-deltaE/Temperature) >= prob_dist(generator_tl))) {
                    // Flip the spin
                    lattice[row][column] *= -1;
                    dE += deltaE;
                    dM += 2.0*lattice[row][column];;
                }
                }
            }


            #pragma omp for 
            for (int row = 1; row<L; row+=2){
                for (int column = 0; column<L; column++){
                // Calculate energy change if this spin is flipped
                double deltaE = calculateDeltaE(lattice, row, column);
                
                // Decide whether to flip the spin based on the Metropolis algorithm
                if (deltaE <= 0 || (exp(-deltaE/Temperature) >= prob_dist(generator_tl))) {
                    // Flip the spin
                    lattice[row][column] *= -1;
                    dE += deltaE;
                    dM += 2.0*lattice[row][column];;
                }
                }
            }


        #pragma omp critical
        {
            energy += dE;
            magnetisation += dM;
        }

            


            #pragma omp barrier
        }
    }


    // cout << energy/(L*L) <<endl;
    // cout << magnetisation/(L*L) << endl;


    double t2 = omp_get_wtime();
    // std::cout << "Elapsed time (seconds): " << (t2 - t1) << std::endl;  // Print elapsed time
    std::cout << (t2 - t1) << std::endl;  // Print elapsed time


    return 0;
}