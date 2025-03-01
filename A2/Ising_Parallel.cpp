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
    int L = 64;
    
    vector<vector<int>> lattice(L, vector<int>(L));

    initializeLattice(lattice);
    
    int sweeps = 10000; // Number of iterations
    double Temperature = 2.5;

    // Get the number of available threads
    int num_threads = 4;
    cout << "Running with " << num_threads << " threads" << endl;

    // Set up separate random number generators for each thread
    vector<mt19937> generators(num_threads);
    vector<uniform_int_distribution<>> position_distribs(num_threads, uniform_int_distribution<>(0, L-1));
    vector<uniform_real_distribution<double>> probability_distribs(num_threads, uniform_real_distribution<double>(0.0, 1.0));
    
    // Initialize each generator with a different seed
    random_device rd;
    for (int i = 0; i < num_threads; i++) {
        generators[i] = mt19937(rd() + i);  // Ensure unique seeds
    }

    // Run simulation with OpenMP
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        mt19937& gen = generators[thread_id];
        uniform_int_distribution<>& pos_dist = position_distribs[thread_id];
        uniform_real_distribution<double>& prob_dist = probability_distribs[thread_id];
        
        // Each thread processes a portion of the total iterations
        #pragma omp for
        for (int i = 0; i < sweeps * L * L; i++) {
            // Choose a random position in the lattice
            int randomRow = pos_dist(gen);
            int randomCol = pos_dist(gen);
            
            // We need critical section when updating the lattice to avoid race conditions
            #pragma omp critical
            {
                // Calculate energy change if this spin is flipped
                double deltaE = calculateDeltaE(lattice, randomRow, randomCol);
                
                // Decide whether to flip the spin based on the Metropolis algorithm
                if (deltaE <= 0 || (exp(-deltaE/Temperature) >= prob_dist(gen))) {
                    // Flip the spin
                    lattice[randomRow][randomCol] *= -1;
                }
            }
        }
    }

    return 0;
}