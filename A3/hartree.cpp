#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>
#include <omp.h>  // Include OpenMP header

using namespace std;

#include <random>
#include "bspline.hpp"
#include "calculateYK.hpp"
#include "functions.hpp"

std::ofstream output_file_1("Hartree_Energy.txt");
std::ofstream output_file_2("Hartree_Wavefunctions.txt");
std::ofstream output_file_3("Hartree_Probability_Density.txt");
std::ofstream output_file_4("Hartree_delta_E.txt");



int main() {

    double r0 = 1.0e-8;
    double rmax = 100.0;
    int num_steps = 10000;
    int k_order = 7;
    int num_splines = 60;
    // 1. Construct Grid
    double dr = (rmax-r0)/num_steps;
    std::vector<double> r(num_steps);
    

    double Z = 3;
    double h = 1;
    double d = 0.2;



    for (int i = 0; i< num_steps; ++i){
        r[i] = r0+dr*i;
    }
    
    std::vector<double> wave_1s(num_steps, 0.0);
    std::vector<double> V_Dir = ykab(0,wave_1s,wave_1s,r);



    // 2. Form B-sline
    auto b_spl = form_Bsplines(r0, rmax, num_steps, k_order, num_splines);
    auto db_spl = form_dBsplines(r0, rmax, num_steps, k_order, num_splines);



    double deltaE_1s = 0;
    double energy1s = 0;
    double newenergy1s = 0;

    for (int iterations=0;iterations < 20; iterations++){

        std::cout << "iter = " << iterations << endl;

        for(int l=0;l<2;++l){

            // 2. Form potential V(r)
            std::vector<double> v(num_steps);
            
            for (int i = 0; i< num_steps; ++i){
                v[i] = -Z/r[i] + l*(l+1)/(2*pow(r[i],2)) + V_Dir[i];
            }
            
            // 3. form H and B
            auto H =  form_H(b_spl, db_spl, v, r0, dr, num_steps);
            auto B =  form_B(b_spl, r0, dr, num_steps);

            auto [EVectors, EValues] = GeneralisedEigenvalue(H, B);

            // After I get EVectors and EValues:
            for(std::size_t n = 0; n < 2; ++n) {
                // Reconstruct the wavefunction for this state

                output_file_1 << iterations << " " << l << " " << n+1 << " " << EValues.at(n) <<"\n";
                std::vector<double> wavefunction(num_steps, 0.0);
                
                // Combine B-splines with eigenvector coefficients to get the wavefunction
                for (int i = 0; i < num_steps; ++i) {
                    for (int j = 0; j < num_splines; ++j) {
                        wavefunction[i] += EVectors(n, j) * b_spl[j][i];
                        
                    }
                }
                
                normalise(wavefunction,num_steps,r,dr);

                for (int i = 0; i < num_steps; ++i) {
                    double radial_prob_density = wavefunction[i] * wavefunction[i] * r[i] * r[i] * 4.0 * M_PI;
                    output_file_3 << iterations << " " << l << " " << n+1 << " " << r[i] << " " << radial_prob_density << "\n";
                    output_file_2 << iterations << " " << l << " " << n+1 << " " << r[i] << " " << wavefunction[i] << "\n";
                }


                
                if (n==0 && l==0){
                    wave_1s = wavefunction;
                    energy1s = newenergy1s;
                    newenergy1s = EValues.at(0);
                    deltaE_1s = newenergy1s-energy1s;
                }

            }
        }


        V_Dir = 2*ykab(0,wave_1s,wave_1s,r);

        if (abs(deltaE_1s) < 1e-6){
            cout << "small energy change" << endl;
            break;
        }

    }

    
    return 0;
}





// v_GR[i] = (Z-1)/r[i] * (h*(exp(r[i]/d)-1))/(1+h*(exp(r[i]/d)-1));


// double VGR_2s = 0;
// double VGR_2p = 0;
// double VEE_2s = 0;
// double VEE_2p = 0;

// double delta_E2s;
// double delta_E2p;

// std::vector<double> v_GR(num_steps,0.0);


// if (l==0 && n==1){

//     for(int i=0; i<num_steps; ++i){
//         VGR_2s += abs(wavefunction[i]) * abs(wavefunction[i]) * v_GR[i] * r[i] * r[i] * 4.0 * M_PI * dr;
//         VEE_2s += abs(wavefunction[i]) * abs(wavefunction[i]) * ykab(0,wave_1s,wave_1s,r)[i] * r[i] * r[i] * 4.0 * M_PI * dr;
//     }
// }
// if (l==1 && n==1){

//     for(int i=0; i<num_steps; ++i){
//         VGR_2p += abs(wavefunction[i]) * abs(wavefunction[i]) * v_GR[i] * r[i] * r[i] * 4.0 * M_PI * dr;
//         VEE_2p += abs(wavefunction[i]) * abs(wavefunction[i]) * ykab(0,wave_1s,wave_1s,r)[i] * r[i] * r[i] * 4.0 * M_PI * dr;
//     }
// }


// delta_E2s = 2*VEE_2s - VGR_2s;
// delta_E2p = 2*VEE_2p - VGR_2p;


// output_file_4 << delta_E2s << endl;
    
// output_file_4 << delta_E2p << endl;