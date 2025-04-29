#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>
#include <vector>

// using namespace std; getting rid of this to avoid name clashes

#include <random>
#include "eigen.hpp"
#include "bspline.hpp"
#include "calculateYK.hpp"
#include "functions.hpp"




int main() {
    
    std::ofstream output_file_1("Greens_Energy.txt");
    std::ofstream output_file_2("Greens_Wavefunctions.txt");
    std::ofstream output_file_3("Greens_Probability_Density.txt");
    std::ofstream output_file_4("Greens_delta_E.txt");

    double r0 = 1.0e-3;
    double rmax = 60.0;
    int num_steps = 4001;
    int k_order = 7;
    int num_splines = 60;
    // 1. Construct Grid
    double dr = (rmax-r0)/num_steps;
    std::vector<double> r(num_steps);
    
    for (int i = 0; i< num_steps; ++i){
        r[i] = r0+dr*i;
    }
    

    // 2. Form B-sline
    auto b_spl = form_Bsplines(r0, rmax, num_steps, k_order, num_splines);
    auto db_spl = form_dBsplines(r0, rmax, num_steps, k_order, num_splines);


    std::vector<double> wave_1s(num_steps, 0.0);
    std::vector<double> wave_2s(num_steps, 0.0);
    std::vector<double> wave_2p(num_steps, 0.0);



    for(int l=0;l<2;++l){

        std::cout << "L = " << l << std::endl;

        // 2. Form potential V(r)
        double Z = 3;
        double h = 1;
        double d = 0.2;
        // double l=0;
        std::vector<double> v(num_steps);
        std::vector<double> v_GR(num_steps, 0.0);

        
        for (int i = 0; i< num_steps; ++i){
            v_GR[i] = (Z-1)/r[i] * (h*(exp(r[i]/d)-1))/(1+h*(exp(r[i]/d)-1));
            v[i] = -Z/r[i] + l*(l+1)/(2*pow(r[i],2)) + v_GR[i];
        }
        

        // 3. form H and B
        auto H =  form_H(b_spl, db_spl, v, r0, dr, num_steps);
        auto B =  form_B(b_spl, r0, dr, num_steps);

        auto [EVectors, EValues] = GeneralisedEigenvalue(H, B);

        for(std::size_t n=l; n<EValues.size(); ++n){
            output_file_1 << l << " " << n+1 << " " << EValues.at(n-l) <<"\n";
        }





        // After I get EVectors and EValues:
        for(std::size_t n = l; n < 2; ++n) {
            // Reconstruct the wavefunction for this state
            std::vector<double> wavefunction(num_steps, 0.0);
            
            make_wavefunction(wavefunction, EVectors, b_spl, num_steps, num_splines, n, l);
            
            normalise(wavefunction,num_steps,r,dr);

            for (int i = 0; i < num_steps; ++i) {
                double radial_prob_density = wavefunction[i] * wavefunction[i];
                output_file_3 << l << " " << n+1 << " " << r[i] << " " << radial_prob_density << "\n";
                output_file_2 << l << " " << n+1 << " " << r[i] << " " << wavefunction[i] << "\n";
            }


            
            if (n==0 && l==0){
                std::cout << "1s wavefunction" << std::endl;
                for (int i = 0; i < num_steps; ++i) {
                    wave_1s[i] = wavefunction[i];
                }
            }
            
            if (l==0 && n==1){
                std::cout << "2s wavefunction" << std::endl;
                double VGR_2s = compute_VGR(wavefunction, v_GR, r, dr);
                double VEE_2s = compute_VEE(wavefunction,wave_1s,r,dr);
                double delta_E2s = VEE_2s - VGR_2s;
                output_file_4 << l << " " <<  delta_E2s << std::endl;
                wave_2s = wavefunction;
            }

            if (l==1 && n==1){
                std::cout << "2p wavefunction" << std::endl;
                double VGR_2p = compute_VGR(wavefunction, v_GR, r, dr);
                double VEE_2p = compute_VEE(wavefunction,wave_1s,r,dr);
                double delta_E2p = VEE_2p - VGR_2p;
                output_file_4 << l << " " << delta_E2p << std::endl;
                wave_2p = wavefunction;
            }

        }

    }



    double time = decay_rate(wave_2s, wave_2p, r, r0, dr, num_steps);

    std::cout << "Time (ns): " << time << std::endl;
    
    return 0;
}