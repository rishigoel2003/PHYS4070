#include <iostream>
int main(){
double r0 = 1.0e-5;
double rmax = 20.0;
int k_spine = 7;
int n_spline = 60;
// order of B-splines
// Initialise the B-spline object

BSpline bspl(k_spine, n_spline, r0, rmax);
// Value of the 1st (index=0) B-spline at r=0
std::cout << bspl.b(0, 0.0) << "\n";
// Value of the 6th (index=5) B-spline at r=1.5 au
std::cout << bspl.b(5, 1.5) << "\n";
// Value of the last (index=N-1) B-spline at r=rmax
std::cout << bspl.b(n_spline- 1, rmax) << "\n";

}