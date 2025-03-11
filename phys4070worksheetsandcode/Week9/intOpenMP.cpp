#include <iostream>
#include <omp.h>
#include <cmath>

int main() {
    double t1 = omp_get_wtime();

    double xmin = 0.0; // Lower bound of integration
    double xmax = 1.0; // Upper bound of integration
    int Npts = 10000;  // Number of points within numerical integral

    double dx = (xmax - xmin) / Npts; // Spacing
    double x, y;

    double integral = 0;
    for (int i=0; i<Npts; i++) {
        for (int j=0; j<Npts; j++) {
            x = i*dx;
            y = j*dx;
            integral += dx*dx / (x*x + y*y + 1);
        }
    }

    std::cout.precision(6);
    std::cout << "Accuracy: " << integral/0.63951 << std::endl;

    double t2 = omp_get_wtime();
    std::cout << "Elapsed time (seconds): " << (t2 - t1) << std::endl;  // Print elapsed time

    return 0;
}