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



void sigma_m_x(Matrix &sig, int i, int N){
    for (int k = 0; k<N; ++k){
        if (k==i){
            sig = sig.tensorProduct(X);
        }
        if (k<i){
            sig = I.tensorProduct(sig);
        }
        if (k>i){
            sig = sig.tensorProduct(I);
        }
    }
}


void sigma_m_n_x(Matrix &sig, int m,int n, int N){
    for (int k = 0; k<N; ++k){
        if (k==m || k==n){
            sig = sig.tensorProduct(X);
        }
        else if (k<min(m,n)){
            sig = I.tensorProduct(sig);
        }
        else if (k>min(m,n)){
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











// Complex vector to represent quantum states
using StateVector = std::vector<std::complex<double>>;

// Function to multiply a matrix by a complex vector: H|ψ⟩
StateVector matrixVectorMultiply(const Matrix& H, const StateVector& psi) {
    size_t size = psi.size();
    StateVector result(size, std::complex<double>(0, 0));
    
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            result[i] += std::complex<double>(H(i, j), 0) * psi[j];
        }
    }
    
    return result;
}






// Compute derivative: dψ/dt = -i·H·ψ (setting ℏ = 1)
StateVector computeDerivative(const Matrix& H, const StateVector& psi) {
    StateVector Hpsi = matrixVectorMultiply(H, psi);
    StateVector dpsi(psi.size());
    
    // Multiply by -i
    for (size_t i = 0; i < psi.size(); i++) {
        dpsi[i] = std::complex<double>(0, -1) * Hpsi[i];
    }
    
    return dpsi;
}



// Perform RK4 step
StateVector rk4Step(const Matrix& H, const StateVector& psi, double dt) {
    size_t size = psi.size();
    
    // Compute k1 = dt * f(t, psi)
    StateVector k1 = computeDerivative(H, psi);
    for (auto& val : k1) val *= dt;
    
    // Compute k2 = dt * f(t + dt/2, psi + k1/2)
    StateVector psi_temp(size);
    for (size_t i = 0; i < size; i++) {
        psi_temp[i] = psi[i] + k1[i] * 0.5;
    }
    StateVector k2 = computeDerivative(H, psi_temp);
    for (auto& val : k2) val *= dt;
    
    // Compute k3 = dt * f(t + dt/2, psi + k2/2)
    for (size_t i = 0; i < size; i++) {
        psi_temp[i] = psi[i] + k2[i] * 0.5;
    }
    StateVector k3 = computeDerivative(H, psi_temp);
    for (auto& val : k3) val *= dt;
    
    // Compute k4 = dt * f(t + dt, psi + k3)
    for (size_t i = 0; i < size; i++) {
        psi_temp[i] = psi[i] + k3[i];
    }
    StateVector k4 = computeDerivative(H, psi_temp);
    for (auto& val : k4) val *= dt;
    
    // Update psi: psi_new = psi + (k1 + 2*k2 + 2*k3 + k4)/6
    StateVector psi_new(size);
    for (size_t i = 0; i < size; i++) {
        psi_new[i] = psi[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    
    return psi_new;
}




// Normalize a state vector
void normalizeState(StateVector& psi) {
    double norm = 0.0;
    for (const auto& val : psi) {
        norm += std::norm(val);
    }
    norm = std::sqrt(norm);
    
    for (auto& val : psi) {
        val /= norm;
    }
}



double dotProduct(const StateVector& a, const StateVector& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("StateVectors must have the same size");
    }
    
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += real(std::conj(a[i]) * b[i]);
    }
    
    return result;
}






double SZ(const StateVector psi,int N){

    double val = 0;
    for (int i = 0; i< N; ++i){
        Matrix sig_i(std::vector<std::vector<double>>{{-1}});
        sigma_m(sig_i,i,N);
        StateVector RHS = matrixVectorMultiply(sig_i,psi);
        val += dotProduct(psi,RHS);
    }

    return val;
}





double SX(const StateVector psi,int N){

    double val = 0;
    for (int i = 0; i< N; ++i){
        Matrix sig_i(std::vector<std::vector<double>>{{-1}});
        sigma_m_x(sig_i,i,N);
        StateVector RHS = matrixVectorMultiply(sig_i,psi);
        val += dotProduct(psi,RHS);
    }

    return val;
}




double CXX(const StateVector psi,int N){

    double val = 0;
    for (int m = 0; m< N; ++m){
        for (int n = 0; n<N; ++n){
            if (n!=m){
                Matrix sig_i(std::vector<std::vector<double>>{{-1}});
                sigma_m_n_x(sig_i,m,n,N);
                StateVector RHS = matrixVectorMultiply(sig_i,psi);
                val += dotProduct(psi,RHS);
            }
        }

    }

    return val;
}







#endif 