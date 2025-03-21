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


    int N = 8;
    int size = std::pow(2, N);  // Size of the matrix (2^N)
    

    int num_steps = 200;
    double dg = 8.0/num_steps;
    
    
    double g = 0;

    // Vectors to store g values and ground state energies
    std::vector<double> g_values(num_steps);
    std::vector<double> ground_state_values(num_steps);
    std::vector<double> theoretical_values(num_steps);



    for (int i=0;i<num_steps; i++){
        g += dg;
        g_values[i] = g;

        theoretical_values[i] = theoreticalGroundStateEnergy(g);

        Matrix I({{1, 0}, {0, 1}}); // Identity matrix
        Matrix Z({{0.5, 0}, {0, -0.5}}); // sigmaZ matrix
        Matrix X({{0, 0.5}, {0.5, 0}}); // sigmaX matrix

        Matrix first = ZTerm(N);
        Matrix second = XTerm(N,g);
        
        Matrix Hamiltonian(size,size);
        Hamiltonian += first;
        Hamiltonian += second;

        auto [eigenvalues, eigenvectors] = Hamiltonian.eigenDecomposition();

        ground_state_values[i] = eigenvalues[0] / N;
    }


    std::vector<double> second_derivatives = calculateSecondDerivative(g_values, ground_state_values);
    std::vector<double> theoretical_second_derivatives = calculateSecondDerivative(g_values, theoretical_values);


    std::ofstream outFile("ground_state.txt");
    outFile << "g\tGround State Energy\tSecond Derivative\tTheory\tTheory Second Derivative" << std::endl;

    for (int i = 0; i < num_steps; i++) {
        outFile << g_values[i] << "\t" 
                << ground_state_values[i] << "\t" 
                << second_derivatives[i] << "\t" 
                << theoretical_values[i] << "\t" 
                << theoretical_second_derivatives[i] << endl;
    }
    
    
    outFile.close();
    std::cout << "Results saved to ground_state.txt" << endl;

    return 0;
}