#ifndef INITIALISE_HPP
#define INITIALISE_HPP

#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>
#include <random>

using namespace std;


//function to initialise lattice with ups and downs
void initializeLattice(vector<vector<int>>& lattice) {
    int L = lattice.size();

    mt19937 gen(69);  // Random number generator (chose a fixed seed)
    uniform_real_distribution<double> dist(0, 1);

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            //want to do a 85 percent bias to one of the spins so we avoid the 50/50 split local minima
            lattice[i][j] = (dist(gen) < 0.85) ? 1 : -1; //ternary operator for if else (if 1 do 1, if 0 do -1)
        }
    }
}



// function to loop through lattic and print the entries of teh matrix
void printLattice(const vector<vector<int>> lattice){
    int L = lattice.size();

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            std::cout << (lattice[i][j] == 1 ? "↑ " : "↓ "); // this is another ternary operator which is shorthand for if else
        }
        std::cout << '\n';
    }

}



// initial energy by taking all neighbours and then dividing by 2
double ComputeInitialEnergy(vector<vector<int>>& lattice) {
    int L = lattice.size();
    double energy = 0; 

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            //sum of all neighbours with periodic boundary conditions
            int neighbors = lattice[(i+1)%L][j] + lattice[i][(j+1)%L] + lattice[(i-1+L)%L][j] + lattice[i][(j-1+L)%L];
            
            energy += -lattice[i][j] * neighbors;
        }
    }
    return energy/2;
}



//initial mag by just taking sum of lattice
double computeInitialMagnetization(const vector<vector<int>>& lattice) {
    int L = lattice.size();
    int sum = 0;
    
    // Sum all spins
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            sum += lattice[i][j];
        }
    }
    
    // Return the mean magnetization per spin
    return sum;
}




#endif // INITIALISE_HPP

