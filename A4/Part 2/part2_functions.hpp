#ifndef PART2_FUNCTIONS_HPP
#define PART2_FUNCTIONS_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#define M_PI 3.14159265358979323846



using namespace std;
using Complex = complex<double>;
using ComplexVector = vector<Complex>;

#include "Matrix.hpp"


Matrix I({{1, 0}, {0, 1}}); // Identity matrix
Matrix Z({{0.5, 0}, {0, -0.5}}); // sigmaZ matrix
Matrix X({{0, 0.5}, {0.5, 0}}); // sigmaX matrix

void sigma_m(Matrix &sig, int i, int N){
    for (int k = 0; k<N; ++k){
        if (k==i){
            sig = sig.tensorProduct(Z);
        }
        if (k<i){
            sig = I.tensorProduct(sig);
        }
        if (k>i){
            sig = sig.tensorProduct(I);
        }
    }
}


void sigma_m_m_plus_1(Matrix &sig, int i, int N){
    for (int k = 0; k<N; ++k){
        if (k==i||k==(i+1)%N){
            sig = sig.tensorProduct(X);
        }
        if (k<i){
            sig = I.tensorProduct(sig);
        }
        if (k>i+1){
            sig = sig.tensorProduct(I);
        }
    }
}




Matrix ZTerm(int N){
    int size = std::pow(2, N);  // Size of the matrix (2^N)
    Matrix first(size,size); // Initialize an empty matrix first

    for (int i=0;i<N;i++){
        Matrix sig_i(std::vector<std::vector<double>>{{-1}});
        sigma_m(sig_i,i,N);
        first += sig_i;
    }
    return first;
}



Matrix XTerm(int N,double g){
    int size = std::pow(2, N);  // Size of the matrix (2^N)
    Matrix interaction_term(size,size); // Initialize an empty matrix first

    for (int i=0;i<N;i++){
        Matrix sig_i_ip1(std::vector<std::vector<double>>{{-g}});
    

        sigma_m_m_plus_1(sig_i_ip1,i,N);
        interaction_term += sig_i_ip1;
    }
    return interaction_term;
}




// Function to calculate second derivative using central difference method
std::vector<double> calculateSecondDerivative(const std::vector<double>& g_values, 
const std::vector<double>& function_values) {
   
    int n = function_values.size();
    std::vector<double> second_derivatives(n, 0.0);

    // For interior points, use central difference
    for (int i = 1; i < n - 1; i++) {
    double h1 = g_values[i] - g_values[i-1];
    double h2 = g_values[i+1] - g_values[i];
    double avg_h = (h1 + h2) / 2.0;

    // Central difference formula for non-uniform grid
    second_derivatives[i] = (function_values[i+1] - 2.0 * function_values[i] + function_values[i-1]) / (avg_h * avg_h);
    }

    // Forward difference for first point (less accurate)
    if (n > 2) {
    double h = g_values[1] - g_values[0];
    second_derivatives[0] = (function_values[2] - 2.0 * function_values[1] + function_values[0]) / (h * h);
    }

    // Backward difference for last point (less accurate)
    if (n > 2) {
    double h = g_values[n-1] - g_values[n-2];
    second_derivatives[n-1] = (function_values[n-1] - 2.0 * function_values[n-2] + function_values[n-3]) / (h * h);
    }

    return second_derivatives;
}






double ellipticIntegralE(double x) {
    const int numIntervals = 1000;  // Number of intervals for integration
    const double upperLimit = M_PI / 2.0;
    const double h = upperLimit / numIntervals;
    
    double sum = 0.0;
    
    // Simpson's rule
    for (int i = 0; i <= numIntervals; i++) {
        double theta = i * h;
        double integrand = std::sqrt(1.0 - x * std::pow(std::sin(theta), 2));
        
        double coefficient = 1.0;
        if (i == 0 || i == numIntervals) {
            coefficient = 1.0;
        } else if (i % 2 == 1) {
            coefficient = 4.0;
        } else {
            coefficient = 2.0;
        }
        
        sum += coefficient * integrand;
    }
    
    return (h / 3.0) * sum;
}


// Function to calculate theoretical ground state energy per site
double theoreticalGroundStateEnergy(double g) {
    double argument = 8.0 * g / std::pow(2.0 + g, 2);
    return -1.0 / (2.0 * M_PI) * (2.0 + g) * ellipticIntegralE(argument);
}







// Function to multiply two matrices
vector<vector<Complex>> multiplyMatrices(const vector<vector<Complex>>& A, const vector<vector<Complex>>& B) {
    size_t n = A.size();
    size_t m = B[0].size();
    size_t k = A[0].size();
    vector<vector<Complex>> result(n, vector<Complex>(m, 0.0));

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            for (size_t l = 0; l < k; l++) {
                result[i][j] += A[i][l] * B[l][j];
            }
        }
    }
    return result;
}

// Function to compute the Hermitian transpose (conjugate transpose) of a matrix
vector<vector<Complex>> hermitianTranspose(const vector<vector<Complex>>& M) {
    size_t n = M.size();
    size_t m = M[0].size();
    vector<vector<Complex>> result(m, vector<Complex>(n));

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            result[j][i] = conj(M[i][j]); // Transpose + Conjugate
        }
    }
    return result;
}

// Function to apply a matrix to a vector
vector<Complex> applyMatrixToVector(const vector<vector<Complex>>& M, const vector<Complex>& v) {
    size_t n = M.size();
    vector<Complex> result(n, 0.0);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < M[i].size(); j++) {
            result[i] += M[i][j] * v[j];
        }
    }
    return result;
}




#endif 