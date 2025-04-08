#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>

using namespace std;

#include <random>
#include "bspline.hpp"
#include "calculateYK.hpp"
#include "functions.hpp"

std::ofstream output_file_1("Hydrogen_Like_Energy.txt");
std::ofstream output_file_2("Hydrogen_Like_Radial.txt");
std::ofstream output_file_3("Hydrogen_Like_Probability_Density.txt");



int main() {

    double r0 = 1.0e-6;
    double rmax = 30.0;
    int num_steps = 2501;
    int k_order = 7;
    int num_splines = 60;
    // 1. Construct Grid
    double dr = (rmax-r0)/num_steps;
    std::vector<double> r(num_steps);
    
    for (int i = 0; i< num_steps; ++i){
        r[i] = r0+dr*i;
    }
    

    // 2. Form B-sline
    auto b_spl =  form_Bsplines(r0, rmax, num_steps, k_order, num_splines);
    auto db_spl =  form_dBsplines(r0, rmax, num_steps, k_order, num_splines);





    for(int l=0;l<2;++l){

        std::cout << "L = " << l << endl;

        // 2. Form potential V(r)
        double Z = 3;
        // double l=0;
        std::vector<double> v(num_steps);
        
        for (int i = 0; i< num_steps; ++i){
            v[i] = -1*Z/r[i] + l*(l+1)/(2*r[i]*r[i]);
        }
        

        // 3. form H and B
        auto H =  form_H(b_spl, db_spl, v, r0, dr, num_steps);
        auto B =  form_B(b_spl, r0, dr, num_steps);

        auto [EVectors, EValues] = GeneralisedEigenvalue(H, B);

        // for(auto en : EValues){
        //     std::cout << en <<"\n";
        // }
        for(std::size_t n=l; n<EValues.size(); ++n){
            output_file_1 << l << " " << n+1 << " " << -EValues.at(n-l) <<"\n";
        }


        // After I get EVectors and EValues:
        for(std::size_t n = l; n < EValues.size(); ++n) {
            // Reconstruct the wavefunction for this state
            std::vector<double> wavefunction(num_steps, 0.0);
            
            // Combine B-splines with eigenvector coefficients to get the wavefunction
            for (int i = 0; i < num_steps; ++i) {
                for (int j = 0; j < num_splines; ++j) {
                    wavefunction[i] += EVectors(n-l, j) * b_spl[j][i];
                }
            }
            
            normalise(wavefunction,num_steps,r,dr);

            for (int i = 0; i < num_steps; ++i) {
                double radial_prob_density = wavefunction[i] * wavefunction[i];
                output_file_3 << l << " " << n+1 << " " << r[i] << " " << radial_prob_density << "\n";
            }
            
            // Calculate expectation value of r
            double r_expectation = calculate_r_expectation(r, dr, wavefunction);
            
            // Output the result
            output_file_2 << l << " " << n+1 << " " << r_expectation << "\n";
        }



    }


    // std::cout << EValues.size() << endl;


    return 0;
}