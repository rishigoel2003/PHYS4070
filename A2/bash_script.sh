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
echo "All tasks completed!"

