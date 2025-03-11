// Note: this won't compile - it's not a full solution
// Just meant to be a scaffold for a solution


//1.  Have the B-splines
//2.  Integrate B-splines
//3.  H and B matrices
//4.  Solve Gen. E-value
#include "bspline.hpp"
#include "matrix.hpp"
#include "eigen.hpp"
#include <vector>

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




int main(){

    double r0 = 1.0e-4;
    double rmax = 30.0;
    int num_steps = 1000;
    int k_order = 7;
    int num_splines = 30;
    // 1. Construct Grid
    double dr = (rmax-r0)/num_steps;
    std::vector<double> r(num_steps);
    
    for (int i = 0; i< num_steps; ++i){
        r[i] = r0+dr*i;
    }
    

    // 2. Form B-sline
    auto b_spl =  form_Bsplines(r0, rmax, num_steps, k_order, num_splines);
    auto db_spl =  form_dBsplines(r0, rmax, num_steps, k_order, num_splines);

    // 2. Form potential V(r)
    double Z = 3;
    double l=0;
    std::vector<double> v(num_steps);
    
    for (int i = 0; i< num_steps; ++i){
        v[i] = -Z/r[i] + l*(l+1)/r[i];
    }
    

    // 3. form H and B
    auto H =  form_H(b_spl, db_spl, v, r0, dr, num_steps);
    auto B =  form_B(b_spl, r0, dr, num_steps);

    auto [EVectors, EValues] = GeneralisedEigenvalue(H, B);

    // for(auto en : EValues){
    //     std::cout << en <<"\n";
    // }
    // for(std::size_t i=0; i<EValues.size(); ++i){
    //     std::cout << EValues.at(i) <<"\n";
    // }

    std::cout << EValues.size();
}

