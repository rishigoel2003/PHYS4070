#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>
#include <omp.h>  // Include OpenMP header

using namespace std;

#include <random>
#include "bspline.hpp"
#include "calculateYK.hpp"
#include "functions.hpp"


int main() {

    int K = 4;     // Order of B-splines (cubic splines)
    int N = 10;    // Number of B-splines
    double r0 = 0.01;  // Start of the range (avoid singularity at r=0)
    double rmax = 10.0; // Reasonable cutoff for hydrogen-like atom

    int num_points = 128;


    // Discretize the interval
    std::vector<double> r = uniform(r0, rmax, num_points);

    BSpline bs(K, N, r0, rmax); // initialise Bsplines

    


    // Compute the B-matrix using the trapezoidal rule
    auto B = compute_B_matrix(bs, N, r0, rmax, num_points);

    // Print the B-matrix
    std::cout << "Overlap matrix B_ij:\n";
    print_matrix(B);

    return 0;
}