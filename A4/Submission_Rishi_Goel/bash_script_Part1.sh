set -e  # Exit on error
#
#
# Compile and run question B1
echo "Compiling and running simulation"
g++ Part1/part1_AB.cpp -o Part1/part1_AB && ./Part1/part1_AB &&
echo "Simulation AB Done"
g++ Part1/part1_C.cpp -o Part1/part1_C && ./Part1/part1_C &&
echo "Simulation C Done"
g++ Part1/part1_D.cpp -o Part1/part1_D && ./Part1/part1_D &&
echo "Simulation D Done"
# #
python3 Part1/plotting_part1_A.py
echo "Plotting question A/B Done"
python3 Part1/plotting_part1_C.py
echo "Plotting question C Done"
python3 Part1/plotting_part1_D.py
echo "Plotting question D Done"
#


echo "All tasks completed!"

