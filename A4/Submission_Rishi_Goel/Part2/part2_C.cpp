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

std::ofstream outFile2("Part2/observables.txt");


int main(){

    vector<double> g_init_vals = {0.0,4.0};
    vector<double> g_vals = {0.0, 1.0, 2.0, 3.0, 4.0};
    
    for (int k = 0; k< g_init_vals.size(); ++k){
        cout << "g_init: " << g_init_vals[k] << endl;

    for (int i = 0; i< g_vals.size(); ++i){

        cout << "g: " << g_vals[i] << endl;
            
        int N = 8;
        int size = std::pow(2, N);  // Size of the matrix (2^N)
        
        //initial g value
        double g_init = g_init_vals[k];
        

        Matrix first = ZTerm(N);
        Matrix second = XTerm(N,g_init);
        
        Matrix Hamiltonian(size,size);
        Hamiltonian += first;
        Hamiltonian += second;

        auto [eigenvalues, eigenvectors] = Hamiltonian.eigenDecomposition();

    

        std::vector<double> psi(N);
        for (size_t i = 0; i < N; i++) {
            psi[i] = eigenvectors(i, 0); // First column
        }


        //evolution g value
        double g = g_vals[i];

        second = XTerm(N,g);
        
        Hamiltonian(size,size);
        Hamiltonian += first;
        Hamiltonian += second;


        double total_time = 5.0;
        double dt = 0.05;
        double t = 0.0;

        StateVector psi_complex(psi.size());
        for (size_t i = 0; i < psi.size(); i++) {
            psi_complex[i] = std::complex<double>(psi[i], 0.0);
        }
        


        while (t < total_time) {
            
            t += dt;
            psi_complex = rk4Step(Hamiltonian, psi_complex, dt);
            
            normalizeState(psi_complex);         

            double s_z = SZ(psi_complex,N);
            double s_x = SX(psi_complex,N);
            double c_xx = CXX(psi_complex,N);

            outFile2 << g_init << " " <<  g << " " << t << " " << s_z << " "<< s_x << " " << c_xx << endl;

        }
        
    }
}

    return 0;
}