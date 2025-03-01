#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
using Complex = complex<double>;
using ComplexVector = vector<Complex>;




// Function to calculate the right-hand side of the Schrödinger equation
ComplexVector calculateRHS(const ComplexVector& psi, double dx, double g) {
    int N = psi.size();
    ComplexVector rhs(N, 0.0);
    
    // Apply the second derivative using the finite difference method
    // and add the nonlinear term
    for (int i = 1; i < N - 1; ++i) {
        // Second derivative: [f(x+dx) - 2f(x) + f(x-dx)]/(dx^2)
        Complex d2psi = (psi[i+1] - 2.0*psi[i] + psi[i-1]) / (dx*dx);
        
        // Nonlinear term: g|ψ|²ψ
        Complex nonlinear = g * norm(psi[i]) * psi[i];
        
        // Combined right-hand side: -d²ψ/dx² + g|ψ|²ψ
        rhs[i] = -Complex(0, 1) * (-d2psi + nonlinear);
    }
    
    // Boundary conditions: ψ(-L/2) = ψ(L/2) = 0
    rhs[0] = 0.0;
    rhs[N-1] = 0.0;
    
    return rhs;
}

// Fourth-order Runge-Kutta method
ComplexVector rungeKutta4(const ComplexVector& psi, double dt, double dx, double g) {
    int N = psi.size();
    ComplexVector k1 = calculateRHS(psi, dx, g);
    
    ComplexVector temp(N);
    for (int i = 0; i < N; ++i) {
        temp[i] = psi[i] + dt / 2.0 * k1[i];
    }
    ComplexVector k2 = calculateRHS(temp, dx, g);
    
    for (int i = 0; i < N; ++i) {
        temp[i] = psi[i] + dt / 2.0 * k2[i];
    }
    ComplexVector k3 = calculateRHS(temp, dx, g);
    
    for (int i = 0; i < N; ++i) {
        temp[i] = psi[i] + dt * k3[i];
    }
    ComplexVector k4 = calculateRHS(temp, dx, g);
    
    ComplexVector psiNew(N);
    for (int i = 0; i < N; ++i) {
        psiNew[i] = psi[i] + dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
    }
    
    return psiNew;
}

// Calculate the norm of the wavefunction (should be conserved)
double calculateNorm(const ComplexVector& psi, double dx) {
    double norm_val = 0.0;
    for (const auto& val : psi) {
        norm_val += norm(val);
    }
    return norm_val * dx;
}

// Function to initialise gaussian pulse
void gaussian(ComplexVector& psi, vector<double> x, double A, double sigma, double x0, double k0, int N, double dx) {
    
    for (int i = 0; i < N; ++i) {
        // Gaussian wavepacket: A * exp(-(x-x0)²/(2*sigma²)) * exp(i*k0*x)
        psi[i] = A * exp(-pow(x[i] - x0, 2) / (2.0 * sigma * sigma)) * 
                 exp(Complex(0, 1) * k0 * x[i]);
    }

    // Normalize the wavefunction
    double initialNorm = calculateNorm(psi, dx);
    for (int i = 0; i < N; ++i) {
        psi[i] /= sqrt(initialNorm);
    }
    
}


vector<double> initialise_grid(int N, int L, double dx){
    // Initialize grid
    vector<double> x(N);
    for (int i = 0; i < N; ++i) {
        x[i] = -L/2 + i * dx;
    }
    return x;

}

#endif