#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "bspline.hpp"
#include "matrix.hpp"
#include "eigen.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <complex>

// b_i(r_j)
std::vector<std::vector<double>> form_Bsplines(double r0, double rmax, int num_steps, int k_order, int num_splines);
std::vector<std::vector<double>> form_dBsplines(double r0, double rmax, int num_steps, int k_order, int num_splines);

double integrate(const std::vector<double> &a, const std::vector<double>&b, double r0, double dr, int num_steps);

Matrix form_H(const std::vector<std::vector<double>>& B_splines, double r0, double dr, int num_steps);

Matrix form_B(const std::vector<std::vector<double>>& B_splines, double r0, double dr, int num_steps);


//*****************************

// r0, rmax, num_steps, dr

std::vector<std::vector<double>> form_Bsplines(double r0, double rmax, int num_steps, int k_order, int num_splines){
   std::vector<std::vector<double>> b(num_splines);
  
   BSpline bspl(k_order, num_splines + 3, r0, rmax);
   double dr = (rmax - r0)/(num_steps-1);
   for(int i=0; i<num_splines; ++i){
       for(int jr=0; jr < num_steps; ++jr){
           double r = r0 + dr*jr;
           b.at(i).push_back(bspl.b(i+2,r));
       }
   }
  
   return b;
}

std::vector<std::vector<double>> form_dBsplines(double r0, double rmax, int num_steps, int k_order, int num_splines){
   std::vector<std::vector<double>> db(num_splines);
  
   BSpline bspl(k_order, num_splines + 3, r0, rmax);
   double dr = (rmax - r0)/(num_steps-1);
   for(int i=0; i<num_splines; ++i){
       for(int jr=0; jr < num_steps; ++jr){
           double r = r0 + dr*jr;
            double temp_db = (bspl.b(i+2,r + dr/2) -  bspl.b(i+2,r - dr/2)) / dr;
            db.at(i).push_back(temp_db);
       }
   }
  
   return db;
}




// f(r)*g(r)
double integrate(const std::vector<double> &a, const std::vector<double>&b, double r0, double dr, int num_steps){
    double integral = 0.0;
    for(int i=0; i<num_steps; ++i){
        integral += a.at(i) * b.at(i) * dr;
    }
    return integral;
}



//f(r)*g(r)*h(r)
double integrate(const std::vector<double> &a, const std::vector<double>&b, const std::vector<double>&c, double r0, double dr, int num_steps){
    double integral = 0.0;
    for(int i=0; i<num_steps; ++i){
        integral += a.at(i) * b.at(i) * c.at(i) * dr;
    }
    return integral;
}


Matrix form_H(const std::vector<std::vector<double>>& B_splines, const std::vector<std::vector<double>>& dB_splines, const std::vector<double> &V, double r0, double dr, int num_steps){
    auto Nspl =  B_splines.size();
    Matrix H(Nspl, Nspl);
    // N_spl * N_spl
    for(int i = 0; i<Nspl; ++i){
        for(int j = 0; j<Nspl; ++j){
            // H_ij = <i|H|j>
            //
            double temp_h1 = 0.5*integrate(dB_splines.at(i), dB_splines.at(j), r0, dr, num_steps);
            double temp_h2 = integrate(B_splines.at(i), V, B_splines.at(j), r0, dr, num_steps);
            H.at(i,j) = temp_h1 + temp_h2;
        }
    }
    return H;
}



Matrix form_B(const std::vector<std::vector<double>>& B_splines, double r0, double dr, int num_steps){
    auto Nspl =  B_splines.size();
    Matrix B(Nspl, Nspl);
    // N_spl * N_spl
    for(int i = 0; i<Nspl; ++i){
        for(int j = 0; j<Nspl; ++j){
            B.at(i,j) =  integrate(B_splines.at(i), B_splines.at(j), r0, dr, num_steps);
        }
    }
    return B;
}



double calculate_r_expectation(const std::vector<double>& r, double dr, 
    const std::vector<double>& wavefunction) {

    double expectation = 0.0;

    for (int i = 0; i < r.size(); ++i) {
        // |ψ(r)|² * r * 4πr²
        double prob_density = wavefunction[i] * wavefunction[i];
        expectation += r[i] * prob_density * r[i] * r[i] * 4.0 * M_PI * dr;
    }

    return expectation;
}







// Overload scalar multiplication (scalar * vector)
std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
}

// Overload scalar multiplication (vector * scalar)
std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    return scalar * vec; // Reuse the previous operator
}


void normalise(std::vector<double> &wavefunction, const double num_steps, const std::vector<double> r, const double dr){
    // Normalize the wavefunction
    double norm = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        norm += abs(wavefunction[i]) * abs(wavefunction[i]) * r[i] * r[i] * 4.0 * M_PI * dr;
    }
    norm = std::sqrt(norm);

    wavefunction = wavefunction * (1/norm);
}




#endif