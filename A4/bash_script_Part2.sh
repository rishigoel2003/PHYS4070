set -e  # Exit on error
#
#
# Compile and run part 2

echo "Compiling and running simulation"
g++ Part2/part2_A.cpp -o Part2/part2_A -llapack -lblas && ./Part2/part2_A &&
echo "Simulation A done"

g++ Part2/part2_B.cpp -o Part2/part2_B -llapack -lblas && ./Part2/part2_B &&
echo "Simulation B done"


g++ Part2/part2_C.cpp -o Part2/part2_C -llapack -lblas && ./Part2/part2_C &&
echo "Simulation C done"


python3 Part2/plotting_part2.py
echo "Plotting Done"


echo "All tasks completed!"

