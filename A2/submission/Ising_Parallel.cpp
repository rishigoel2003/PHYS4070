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

std::ofstream output_file("IsingParallel.txt", std::ios::app);


// Main function
int main() {
    double t1 = omp_get_wtime();

    int L = 128;
    int N = pow(L,2);
    double B=0;

    vector<vector<int>> lattice(L, vector<int>(L));

    initializeLattice(lattice);
    double energy = ComputeInitialEnergy(lattice);
    double magnetisation = computeInitialMagnetization(lattice);
    
    int sweeps = 2000; // Number of iterations
    int numTemperatures = 20;
    double minTemp = 0.3; // Avoiding T=0 for stability
    double maxTemp = 5;
    int num_threads = 0;

    for (int tempIndex = 0; tempIndex < numTemperatures; tempIndex++) {
        // Calculate current temperature
        double Temperature = minTemp + tempIndex * (maxTemp - minTemp) / (numTemperatures - 1);
    

    // Run simulation with OpenMP
    #pragma omp parallel default(none) shared(energy, magnetisation, Temperature, sweeps, lattice, L, N,B, num_threads, output_file, std::cout)
    {

        int thread_id = omp_get_thread_num();

        if (thread_id+1 > num_threads) {num_threads = thread_id + 1;}

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
                double deltaE = calculateDeltaE(lattice, row, column,B);
                
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
                double deltaE = calculateDeltaE(lattice, row, column,B);
                
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


            #pragma omp master
            output_file << i << " "<< num_threads << " " << Temperature << " " << energy/N << " " << magnetisation/N << "\n";
            

        }
    }
    }


    // cout << energy/(L*L) <<endl;
    // cout << magnetisation/(L*L) << endl;


    double t2 = omp_get_wtime();
    // std::cout << "Elapsed time (seconds): " << (t2 - t1) << std::endl;  // Print elapsed time
    std::cout << (t2 - t1) << std::endl;  // Print elapsed time

    output_file.close();


    return 0;
}