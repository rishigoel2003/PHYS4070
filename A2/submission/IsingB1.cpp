//packages
#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>

// making the output file which we save data to
std::ofstream output_file("IsingB1.txt");

// i know im not meant to do this but it makes it easier for me to code cout
using namespace std;

#include <random>

//helper functions and initialising matrix
#include "initialise.hpp"
#include "functions.hpp"



// Main function
int main() {
    //params
    int L = 16;
    double B = 0;
    int N = pow(L,2);

    //making and initialising lattice with energy and magnetisation and printing it out
    vector<vector<int>> lattice(L, vector<int>(L));

    initializeLattice(lattice);

    double Temperature = 2;


    double energy = ComputeInitialEnergy(lattice);
    double magnetisation = computeInitialMagnetization(lattice);
    cout << "Initial energy: " << energy << " Initial Magnetisation: " << magnetisation << endl; 

    // can print out if you want
    // printLattice(lattice);

    //making mercer twin random generator
    random_device rd;
    mt19937 gen(rd());
    
    uniform_int_distribution<> distrib(0, L-1);
    uniform_real_distribution<double> dist(0.0, 1.0);



    cout << "Running simulation at T = " << Temperature << endl;
    int sweeps = 2000; // Number of iterations
    



    // Run simulation
    for (int i = 0; i < sweeps*N; i++) {


        // choose a random position in the lattice
        int randomRow = distrib(gen);
        int randomCol = distrib(gen);
        
        // calculate energy change if this spin is flipped
        double deltaE = calculateDeltaE(lattice, randomRow, randomCol, B);
        
        // decide whether to flip the spin based on the Metropolis algorithm (energy less than 0 or with probability e^-dE/T)
        if (deltaE <= 0 || (exp(-deltaE/Temperature) >= dist(gen))) {
            // Flip the spin
            lattice[randomRow][randomCol] *= -1;
            // Update energy
            energy += deltaE;
            //update magnetisation (this is because the magnetisation is independent of the neighbours its just the sum of spins)
            magnetisation += 2.0*lattice[randomRow][randomCol];
        }
        

        if (i % (L*L) == 0) {
            // Take measurements here (i.e. at the end of each sweep)
            output_file << i/(L*L) << " " << Temperature << " " << energy/N << " " << magnetisation/N << "\n";
        }
    
    }
    
    //output details so i know its done
    cout << "B1 algorithm finished\n";

    printLattice(lattice);
    
    cout << "Final energy: " << energy << " Final Magnetisation: " << magnetisation << endl; 
    
    output_file.close();

    return 0;
}





