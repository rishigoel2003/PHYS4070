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

// const std::complex<double> I(0, 1); // Imaginary unit

std::ofstream outFile("part2/psi.txt");
std::ofstream outFile2("part2/observables.txt");


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

    g = 4;

    second = XTerm(N,g);
    
    Hamiltonian(size,size);
    Hamiltonian += first;
    Hamiltonian += second;

    // auto [eigenvalues2, eigenvectors2] = Hamiltonian.eigenDecomposition();

    //need to extend my stuff to complex numbers or use RK4 to evolve wavefunction

    double total_time = 10;
    double dt = 0.1;
    double t = 0.0;

    StateVector psi_complex(psi.size());
    for (size_t i = 0; i < psi.size(); i++) {
        psi_complex[i] = std::complex<double>(psi[i], 0.0);
        outFile <<  t << " " << i << " " << psi_complex[i] << std::endl;
    }
    


    while (t < total_time) {
        
        t += dt;
        psi_complex = rk4Step(Hamiltonian, psi_complex, dt);
        
        normalizeState(psi_complex);
        
        // Output results
        for (size_t i = 0; i < psi_complex.size(); i++) {
            outFile <<  t << " " << i << " " << psi_complex[i] << std::endl;
        }                

        double s_z = SZ(psi_complex,N);
        double s_x = SX(psi_complex,N);
        double c_xx = CXX(psi_complex,N);

        outFile2 << t << " " << s_z << " "<< s_x << " " << c_xx << endl;

    }
    

    return 0;
}