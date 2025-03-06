#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include "bspline.hpp"

using namespace Range;

double integrate_trapezoidal(const std::vector<double>& r, const std::vector<double>& f_values) {
    double integral = 0.0;
    int N = r.size();
    for (int i = 0; i < N - 1; i++) {
        double h = r[i + 1] - r[i];
        integral += 0.5 * h * (f_values[i] + f_values[i + 1]);
    }
    return integral;
}



struct BSplineData {
    std::vector<std::vector<double>> B_vals;   // B-spline values
    std::vector<std::vector<double>> B_derivs; // B-spline derivatives
    std::vector<std::vector<double>> B_2nd_derivs; // B-spline 2nd derivatives
    std::vector<double> r_grid;                // Discretized r values
};



// Compute B-spline values and derivatives on a grid
BSplineData precompute_bspline_data(const BSpline& bs, int N, double r_min, double r_max, int num_points) {
    BSplineData data;
    data.r_grid.resize(num_points);
    data.B_vals.resize(N, std::vector<double>(num_points, 0.0));
    data.B_derivs.resize(N, std::vector<double>(num_points-1, 0.0));
    data.B_2nd_derivs.resize(N, std::vector<double>(num_points-2, 0.0));


    data.r_grid = uniform(r_min,r_max,num_points);
    
    for (int k = 0; k < num_points; k++) {
        for (int i = 0; i < N; i++) {
            data.B_vals[i][k] = bs(i, data.r_grid[k]);
            
        }
    }

    for (int k = 0; k < num_points-1; k++) {
        double dr = data.r_grid[k+1] - data.r_grid[k];
        for (int i = 0; i < N; i++) {
            data.B_derivs[i][k] = (data.B_vals[i][data.r_grid[k] + dr/2] - data.B_vals[i][data.r_grid[k] - dr/2]) / (dr); // Central difference approx
        }
    }

    for (int k = 0; k < num_points-2; k++) {
        double dr = data.r_grid[k+1] - data.r_grid[k];
        for (int i = 0; i < N; i++) {
            data.B_2nd_derivs[i][k] = (data.B_derivs[i][data.r_grid[k] + dr/2] - data.B_derivs[i][data.r_grid[k] - dr/2]) / (dr); // Central difference approx
        }
    }


    return data;
}






// std::vector<std::vector<double>> compute_B_matrix(const BSpline& bs, int N, double r_min, double r_max, int num_points) {
//     std::vector<std::vector<double>> B(N, std::vector<double>(N, 0.0));



//     // Compute B_ij
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j <= i; j++) { // Use symmetry
//             std::vector<double> f_values(num_points);
//             for (int k = 0; k < num_points; k++) {
//                 f_values[k] = bs(i, r[k]) * bs(j, r[k]);
//             }
//             B[i][j] = integrate_trapezoidal(r, f_values);
//             B[j][i] = B[i][j]; // Symmetric matrix
//         }
//     }
//     return B;
// }


// Print matrix
void print_matrix(const std::vector<std::vector<double>>& M) {
    for (const auto& row : M) {
        for (double val : row) {
            printf("%8.5f ", val);
        }
        std::cout << "\n";
    }
}



#endif