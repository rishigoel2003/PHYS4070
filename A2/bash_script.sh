set -e  # Exit on error
#
#
# Compile and run question B1
echo "Compiling and running ising model"
g++ Ising_Basic.cpp -o Ising && ./Ising &&
#
echo "Compiling and running ising model with scaling"
g++ Ising_Scaling.cpp -o Ising_scaling && ./Ising_scaling &&
# Plotting: Shows 3 sets of 4 plots, first set corresponds to B1, 2nd set corresponds to B2 and 3rd corresponds to B2 with kick
# Need to close each set of 4 for the next set to show up
# Not all plots are shown in the assignment report but they are all in the code for completeness
# Plots do not save automatically, have to manually save if you want
echo "Generating plots for first question..."
python3 plotting.py
#
#
echo "Generating plots for second question..."
python3 plotting_scaling.py
#
echo "All tasks completed!"

