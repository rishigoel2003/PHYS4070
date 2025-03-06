#ifndef PART2_FUNCTIONS_HPP
#define PART2_FUNCTIONS_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

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
        if (k==i||k==i+1){
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

    for (int i=0;i<N-1;i++){
        Matrix sig_i_ip1(std::vector<std::vector<double>>{{-g}});
        sigma_m_m_plus_1(sig_i_ip1,i,N);
        interaction_term += sig_i_ip1;
    }
    return interaction_term;
}









#endif 