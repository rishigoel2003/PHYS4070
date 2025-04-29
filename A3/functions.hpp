#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "bspline.hpp"
#include "matrix.hpp"
#include "eigen.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <complex>
#include <algorithm> 


// b_i(r_j)
std::vector<std::vector<double>> form_Bsplines(double r0, double rmax, int num_steps, int k_order, int num_splines);
std::vector<std::vector<double>> form_dBsplines(double r0, double rmax, int num_steps, int k_order, int num_splines);

double integrate(const std::vector<double> &a, const std::vector<double>&b, double r0, double dr, int num_steps);

Matrix form_H(const std::vector<std::vector<double>>& B_splines, double r0, double dr, int num_steps);

Matrix form_B(const std::vector<std::vector<double>>& B_splines, double r0, double dr, int num_steps);



// Overload vector multiplication (vector * vector)
inline std::vector<double> operator*(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] * vec2[i];
    }
    return result;
}




// Overload scalar multiplication (scalar * vector)
inline std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
}

// Overload scalar multiplication (vector * scalar)
inline std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    return scalar * vec; // Reuse the previous operator
}



//*****************************

// r0, rmax, num_steps, dr


// making B Spline matrix
inline std::vector<std::vector<double>> form_Bsplines(double r0, double rmax, int num_steps, int k_order, int num_splines){
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




// making the derivative of B Spline matrix
inline std::vector<std::vector<double>> form_dBsplines(double r0, double rmax, int num_steps, int k_order, int num_splines){
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




// Simpson's rule integration for f(r)*g(r)
inline double integrate(const std::vector<double> &a, const std::vector<double>&b, double r0, double dr, int num_steps) {
    // Ensure odd number of steps for Simpson's rule
    if (num_steps % 2 == 0) {
        throw std::invalid_argument("Number of steps must be odd for Simpson's rule");
    }

    double integral = 0.0;
    
    // First and last points have weight 1
    integral += a.at(0) * b.at(0);
    integral += a.at(num_steps-1) * b.at(num_steps-1);

    // Even-indexed points (except first and last) have weight 4
    for (int i = 1; i < num_steps-1; i += 2) {
        integral += 4 * a.at(i) * b.at(i);
    }

    // Odd-indexed points have weight 2
    for (int i = 2; i < num_steps-1; i += 2) {
        integral += 2 * a.at(i) * b.at(i);
    }

    // Multiply by (dr/3)
    return integral * dr / 3.0;
}




// Simpson's rule integration for f(r)*g(r)*h(r)
inline double integrate(const std::vector<double> &a, const std::vector<double>&b, const std::vector<double>&c, double r0, double dr, int num_steps) {
    // Ensure odd number of steps for Simpson's rule
    if (num_steps % 2 == 0) {
        throw std::invalid_argument("Number of steps must be odd for Simpson's rule");
    }

    double integral = 0.0;
    
    // First and last points have weight 1
    integral += a.at(0) * b.at(0) * c.at(0);
    integral += a.at(num_steps-1) * b.at(num_steps-1) * c.at(num_steps-1);

    // Even-indexed points (except first and last) have weight 4
    for (int i = 1; i < num_steps-1; i += 2) {
        integral += 4 * a.at(i) * b.at(i) * c.at(i);
    }

    // Odd-indexed points have weight 2
    for (int i = 2; i < num_steps-1; i += 2) {
        integral += 2 * a.at(i) * b.at(i) * c.at(i);
    }

    // Multiply by (dr/3)
    return integral * dr / 3.0;
}




// Making the Hamiltonian matrix by integrating the B splines and the potential for each term
inline Matrix form_H(const std::vector<std::vector<double>>& B_splines, const std::vector<std::vector<double>>& dB_splines, const std::vector<double> &V, double r0, double dr, int num_steps){
    auto Nspl =  B_splines.size();
    Matrix H(Nspl, Nspl);
    // finding each term of H
    for(int i = 0; i<Nspl; ++i){
        for(int j = 0; j<Nspl; ++j){
            
            //first term of hamiltonian after integration by parts
            double h1 = 0.5*integrate(dB_splines.at(i), dB_splines.at(j), r0, dr, num_steps);
            //2nd term with potential
            double h2 = integrate(B_splines.at(i), V, B_splines.at(j), r0, dr, num_steps);
            H.at(i,j) = h1 + h2;
        }
    }
    return H;
}



//Modified Hamiltonian matrix for Hartree-Fock method
inline Matrix form_H_Hartree_Fock(const std::vector<std::vector<double>>& B_splines, const std::vector<std::vector<double>>& dB_splines, const std::vector<double> &V, double r0, double dr, int num_steps
    , std::vector<double> wave_1s, std::vector<double> r,int l){
    auto Nspl =  B_splines.size();
    Matrix H(Nspl, Nspl);
    double Lambda = 0;
    if (l==0){Lambda = 1/2.0;}
    if (l==1){Lambda = 1/6.0;}

    for(int i = 0; i<Nspl; ++i){
        for(int j = 0; j<Nspl; ++j){
            
            
            std::vector<double> v_exch = -2*Lambda*ykab(l,wave_1s,B_splines.at(j),r)*wave_1s;

            double h1 = 0.5*integrate(dB_splines.at(i), dB_splines.at(j), r0, dr, num_steps);
            double h2 = integrate(B_splines.at(i), V, B_splines.at(j), r0, dr, num_steps);

            //exchange term
            double h3 = integrate(B_splines.at(i), v_exch, r0, dr, num_steps);
            H.at(i,j) = h1 + h2 + h3;
        }
    }
    return H;
}







//making the B matrix by integrating the B splines (not dirac deltas as they arent normalised and orthogonal)
inline Matrix form_B(const std::vector<std::vector<double>>& B_splines, double r0, double dr, int num_steps){
    auto Nspl =  B_splines.size();
    Matrix B(Nspl, Nspl);
    for(int i = 0; i<Nspl; ++i){
        for(int j = 0; j<Nspl; ++j){
            B.at(i,j) =  integrate(B_splines.at(i), B_splines.at(j), r0, dr, num_steps);
        }
    }
    return B;
}



inline double calculate_r_expectation(std::vector<double>& r, double dr, 
    const std::vector<double>& wavefunction) {

    double expectation = 0.0;
    int num_steps = size(r);
    double r0 = r[0];

    //order doesnt matter for numerical integration
    expectation = integrate(wavefunction,wavefunction,r,r0,dr,num_steps);

    return expectation;
}



inline void normalise(std::vector<double> &wavefunction, const double num_steps, const std::vector<double> r, const double dr){
    // Normalize the wavefunction

    double r0 = r[0];


    // Take absolute value of wavefunction element-wise just in case
    std::vector<double> abs_wavefunction(wavefunction.size());
    std::transform(wavefunction.begin(), wavefunction.end(), abs_wavefunction.begin(), [](double val) {
        return std::abs(val);
    });

    double norm = std::sqrt(integrate(abs_wavefunction,abs_wavefunction,r0,dr,num_steps));

    //normalise
    wavefunction = wavefunction * (1/norm);
}




//computing the perturbation for VGR
inline double compute_VGR(const std::vector<double>& wavefunction, const std::vector<double>& v_GR, const std::vector<double>& r, double dr) {

    int num_steps = size(r);
    double r0 = r[0];

    // Take absolute value of wavefunction element-wise
    std::vector<double> abs_wavefunction(wavefunction.size());
    std::transform(wavefunction.begin(), wavefunction.end(), abs_wavefunction.begin(), [](double val) {
        return std::abs(val);
    });

    double VGR = integrate(abs_wavefunction,abs_wavefunction,v_GR,r0,dr,num_steps);
    return VGR;
}



//computing the perturbation for VEE part of the first order perturbation 
inline double compute_VEE(const std::vector<double>& wavefunction, const std::vector<double>& wave_1s, const std::vector<double>& r, double dr) {


    int num_steps = size(r);
    double r0 = r[0];

    // Take absolute value of wavefunction element-wise
    std::vector<double> abs_wavefunction(wavefunction.size());
    std::transform(wavefunction.begin(), wavefunction.end(), abs_wavefunction.begin(), [](double val) {
        return std::abs(val);
    });

    double VEE = 2*integrate(abs_wavefunction,abs_wavefunction,ykab(0,wave_1s,wave_1s,r),r0,dr,num_steps);
    
    return VEE;
}



//this is so that we can read the wavefunction from the file and get the last iteration of the wavefunction and use in HF
inline std::vector<double> read_wave1s_from_file(const std::string& filename) {
    std::ifstream input_file(filename);

    std::vector<double> wave_1s;
    std::string line;
    int current_iterations = 0, max_iterations = 0;
    int current_l, current_n;
    double r, wavefunction;

    // First pass: find the maximum iteration number
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        if (iss >> current_iterations >> current_l >> current_n >> r >> wavefunction) {
            max_iterations = std::max(max_iterations, current_iterations);
        }
    }

    // Reset file stream to beginning
    input_file.clear();
    input_file.seekg(0, std::ios::beg);

    // Second pass: collect data from the final iteration
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        
        if (!(iss >> current_iterations >> current_l >> current_n >> r >> wavefunction)) {
            throw std::runtime_error("Error parsing line: " + line);
        }

        // Check if this line is from the final iteration and matches desired l and n
        if (current_iterations == max_iterations && current_l == 0 && current_n == 1) {
            wave_1s.push_back(wavefunction);
        }
    }



    return wave_1s;
}


inline std::vector<double> make_wavefunction(std::vector<double>& wavefunction, const Matrix& EVectors, const std::vector<std::vector<double>>& b_spl, int num_steps, int num_splines, int n, int l) {
    // Reconstruct the wavefunction for this state
    wavefunction.resize(num_steps, 0.0);
// Combine B-splines with eigenvector coefficients to get the wavefunction
for (int i = 0; i < num_steps; ++i) {
    for (int j = 0; j < num_splines; ++j) {
        wavefunction[i] += EVectors(n-l, j) * b_spl[j][i];
    }
}
return wavefunction;
}


inline double decay_rate(std::vector<double>& wave_2s, std::vector<double>& wave_2p, const std::vector<double>& r, double r0, double dr, int num_steps) {
    // Calculate the decay rate of the 2s to 2p transition
    // using the formula: decay_rate = 2 * R_ab^2 * omega_ab^3 / 3
    double R_ab = integrate(wave_2s, wave_2p, r, r0, dr, num_steps);
    double omega_ab = 0.06791;
    double unit_change = 1.071e10;
    double decay_rate = 2 * std::pow(R_ab, 2) * std::pow(omega_ab, 3) / 3 * unit_change;

    double time = std::round((1 / decay_rate * 1e9) * 1e4) / 1e4; // rounding to 4 decimal places

    return time;
}





#endif