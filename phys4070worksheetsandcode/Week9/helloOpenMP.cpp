// To compile, add the fopenmp flag:
// g++ helloopenmp.cpp -o helloopenmp -fopenmp
// Before running the code, choose the number of threads by type into the terminal:
// export OMP_NUM_THREADS=4

// Run it a few times. You should see that the output is non-deterministic

#include <iostream>
#include <omp.h>

int main()
{
	// Declare "nthreads" before the parallel region: this variable will be shared by all threads
	int nthreads;

	nthreads = omp_get_num_threads(); // This in-built function returns the total number of threads
	std::cout << "Number of threads = " << nthreads << std::endl; // Print the total number of threads

	// This next line creates a parallel region: a team of threads is created, and each thread runs this next section of code in parallel
	// Important: it's best practice to use "default(none)", and declare all private variables within the parallel region
	// Inside "shared()", list all the variables you want to be shared across all threads (cout is a variable, so must be included if you want output)
	#pragma omp parallel default(none), shared(nthreads, std::cout)
	{ // Start of parallel region: everything inside these braces will be run by multiple threads
		
		// This variable is declared inside the parallel region, so it will be private to each thread
		int ntid = omp_get_thread_num(); // This in-built function returns this thread's index (0, 1, 2, 3...)
		std::cout << "Hello world from thread " << ntid << std::endl; // Print thread ID

		// Only "master" thread (thread 0) does this next part
		if (ntid == 0){
			
			nthreads = omp_get_num_threads(); // Get the total number of threads again
			std::cout << "Number of threads = " << nthreads << std::endl; // Print
		}

	}  // End of parallel region: all threads join master thread and terminate

	return 0;
}