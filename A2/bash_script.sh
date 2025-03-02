set -e  # Exit on error
#
#
# Compile and run question B1
echo "Compiling and running ising model"
g++ Ising_Basic.cpp -o Ising && ./Ising &&
#
echo "Compiling and running ising model with scaling"
g++ Ising_Scaling.cpp -o Ising_scaling && ./Ising_scaling &&
#
#
echo "Compiling and running ising model for power laws"
g++ Ising_Power_Laws.cpp -o Ising_Power_Laws && ./Ising_Power_Laws &&
#
echo "Generating plots for first question..."
python3 plotting.py
#
#
echo "Generating plots for second question..."
python3 plotting_scaling.py
#
echo "Generating plots for third question..."
python3 plotting_power.py
#
#
#
#
#
# Output file for timing results
OUTPUT_FILE="timing_results.txt"
# Clear previous results
echo "Threads,Time" > $OUTPUT_FILE
# Compile the program
g++ Ising_Parallel.cpp -o parallel -fopenmp
# Run with different thread counts
for threads in 1 2 4 8 16
do
  echo "Running with $threads threads..."
  export OMP_NUM_THREADS=$threads
  # Run the program and capture the time output directly
  time_value=$(./parallel)
  # Append to results file
  echo "$threads,$time_value" >> $OUTPUT_FILE
#   echo "Completed run with $threads threads: $time_value seconds"
done
#plotting for the parallelisation question using OpenMP
echo "Generating plots for OpenMP question..."
python3 plot_times.py
#
#
#
echo "All tasks completed!"