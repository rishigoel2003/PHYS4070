set -e  # Exit on error
# #
# #
#making a plots folder for the plots to go into
if [ ! -d "Plots" ]; then
  mkdir -p Plots
  echo "Created Plots directory"
fi
# Compile and run question B1
echo "Compiling and running ising model"
g++ IsingB1.cpp -o B1 && ./B1 &&
#
echo "Compiling and running ising model with scaling"
g++ Ising_Scaling.cpp -o Ising_scaling && ./Ising_scaling &&
#
echo "Compiling and running ising model with applied mag field"
g++ Ising_AppliedB.cpp -o Ising_AppliedB && ./Ising_AppliedB &&
#
echo "Compiling and running ising model for power laws"
g++ Ising_Power_Laws.cpp -o Ising_Power_Laws && ./Ising_Power_Laws &&
#
#
#
echo "Generating plots for first question..."
python3 B1.py
#
#
echo "Generating plots for second question..."
python3 plotting_scaling.py
#
echo "Generating plots for third question..."
python3 plotting_AppliedB.py
#
echo "Generating plots for fourth question..."
python3 plotting_power.py
#
#
#
#
#delete this file so we can remake it and add data to it fresh
rm -f IsingParallel.txt
# Output file for timing results
OUTPUT_FILE="timing_results.txt"
# put in header to file
echo "Threads,Time" > $OUTPUT_FILE
# Compile the program with openmp
g++ Ising_Parallel.cpp -o parallel -fopenmp
# Run with different thread counts
for threads in 1 2 4 8 16
do
  echo "Running with $threads threads..."
  export OMP_NUM_THREADS=$threads
  # Run the program and capture the time output directly as we are only outputting a single value which is the time
  time_value=$(./parallel)
  # Append to results file
  echo "$threads,$time_value" >> $OUTPUT_FILE
  echo "Completed run with $threads threads: $time_value seconds"
done
#plotting for the parallelisation question using OpenMP
echo "Generating plots for OpenMP question..."
python3 plot_times.py
#
#
#
echo "All tasks completed!"