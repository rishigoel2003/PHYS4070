#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <chrono>

using namespace std;
using Complex = complex<double>;
using ComplexVector = vector<Complex>;

#include "Matrix.hpp"
#include "part2_functions.hpp"



int main() {
    // Open file for writing results
    std::ofstream outFile("execution_times.txt");
    
    // Write header to file
    outFile << "N\tMatrix Size\tExecution Time (ms)" << std::endl;
    
    for (int N = 2; N < 9; N++) {
        auto start = std::chrono::high_resolution_clock::now();

        int size = std::pow(2, N);  // Size of the matrix (2^N)
        double g = 1.0;
        
        // Start timing
        
        Matrix I({{1, 0}, {0, 1}}); // Identity matrix
        Matrix Z({{0.5, 0}, {0, -0.5}}); // sigmaZ matrix
        Matrix X({{0, 0.5}, {0.5, 0}}); // sigmaX matrix
        
        Matrix first = ZTerm(N);
        Matrix second = XTerm(N, g);
        
        Matrix Hamiltonian(size, size);
        Hamiltonian += first;
        Hamiltonian += second;
        
        auto [eigenvalues, eigenvectors] = Hamiltonian.eigenDecomposition();
        
        // End timing
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        // Write results to file
        outFile << N << "\t" << size << "\t" << duration << std::endl;
        
        // Also output to console for monitoring progress
        std::cout << "N = " << N << ", Matrix Size = " << size 
                  << ", Time = " << duration << " ms" << std::endl;
    }
    
    outFile.close();
    std::cout << "Results saved to execution_times.txt" << std::endl;
    
    return 0;
}