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



int main(){


    int N = 4;
    int size = std::pow(2, N);  // Size of the matrix (2^N)
    double g = 1.0;

    Matrix I({{1, 0}, {0, 1}}); // Identity matrix
    Matrix Z({{0.5, 0}, {0, -0.5}}); // sigmaZ matrix
    Matrix X({{0, 0.5}, {0.5, 0}}); // sigmaX matrix

    Matrix first = ZTerm(N);
    Matrix second = XTerm(N,g);
    
    Matrix Hamiltonian(size,size);
    Hamiltonian += first;
    Hamiltonian += second;

    // Hamiltonian.print();

    auto [eigenvalues, eigenvectors] = Hamiltonian.eigenDecomposition();

    // Print results
    std::cout << "Eigenvalues:\n";
    for (auto val : eigenvalues) {
        std::cout << val << " ";
    }

    
    std::cout << "\n\nEigenvectors (as columns):\n";
    eigenvectors.print();


    return 0;
}