set -e  # Exit on error
#
#
# Compile and run question B1
echo "Compiling and running simulation"
g++ part1_AB.cpp -o part1_AB && ./part1_AB &&
g++ part1_C.cpp -o part1_C && ./part1_C &&
g++ part1_D.cpp -o part1_D && ./part1_D &&
#
python3 plotting_part1.py
#
echo "All tasks completed!"

