// To compile, add the fopenmp flag:
// g++ sumopenmp.cpp -o sumopenmp -fopenmp
// Before running the code, choose the number of threads by type into the terminal:
// export OMP_NUM_THREADS=4

#include <iostream>
#include <omp.h>

int main()
{
    double t1 = omp_get_wtime(); // Start time
	int N = 10000; // Number of iterations
    //double A = 0; // Variable to sum over

    // Start parallel region
    #pragma omp parallel default(none) shared(N, std::cout)//, A)
    {   
        int ntid = omp_get_thread_num(); // Private variable: thread index
        int count = 0; // Private variable: iteration counter for this thread

        // Parallelise the following for loop
        #pragma omp for
    	for (int i=0; i<N; i++) {

            //A += 1;

            // Add to this thread's counter
            count++;
    	}

        // Print output from each thread. Critical directive makes sure only one thread runs this line at a time (just to make the output clean)
        #pragma omp critical
        {
            std::cout << "Thread " << ntid << " ran a total of " << count << " iterations (" << 100.0*(double)count/(double)N << "%)" << std::endl;
        }
    }

    double t2 = omp_get_wtime(); // End time
    std::cout << "Elapsed time (seconds): " << (t2 - t1) << std::endl;  // Print elapsed time
    //cout << A << endl;

	return 0;
}