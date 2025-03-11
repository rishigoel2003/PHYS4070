#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>
#include <random>

using namespace std;


double calculateDeltaE(const vector<vector<int>>& lattice, int row, int col,double B) {
    int L = lattice.size();
    int spin = lattice[row][col];
    
    // Get the spin's neighbors (with periodic boundary conditions)
    int top = lattice[((row - 1) + L) % L][col];
    int bottom = lattice[((row + 1) + L) % L][col];
    int left = lattice[row][((col - 1) + L) % L];
    int right = lattice[row][((col + 1) + L) % L];
    
    // Calculate current interaction energy
    int sum_neighbors = top + bottom + left + right;
    
    // Energy change if we flip this spin
    return 2 * spin * (sum_neighbors + B);
}








#endif // FUNCTIONS_HPP

