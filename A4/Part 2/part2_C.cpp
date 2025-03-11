#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdexcept>



using namespace std;
using Complex = complex<double>;
using ComplexVector = vector<Complex>;


#include "Matrix.hpp"
#include "part2_functions.hpp"

const std::complex<double> I(0, 1); // Imaginary unit



int main(){


    int N = 8;
    int size = std::pow(2, N);  // Size of the matrix (2^N)
    
    double g = 0;
    
    Matrix I({{1, 0}, {0, 1}}); // Identity matrix
    Matrix Z({{0.5, 0}, {0, -0.5}}); // sigmaZ matrix
    Matrix X({{0, 0.5}, {0.5, 0}}); // sigmaX matrix

    Matrix first = ZTerm(N);
    Matrix second = XTerm(N,g);
    
    Matrix Hamiltonian(size,size);
    Hamiltonian += first;
    Hamiltonian += second;

    auto [eigenvalues, eigenvectors] = Hamiltonian.eigenDecomposition();

    std::vector<double> psi(N);
    for (size_t i = 0; i < N; i++) {
        psi[i] = eigenvectors(i, 0); // First column
    }

    double g = 4;

    Matrix second = XTerm(N,g);
    
    Matrix Hamiltonian(size,size);
    Hamiltonian += first;
    Hamiltonian += second;

    auto [eigenvalues, eigenvectors] = Hamiltonian.eigenDecomposition();

    //need to extend my stuff to complex numbers



    return 0;
}