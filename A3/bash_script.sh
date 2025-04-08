set -e  # Exit on error

echo "Making Plots Directory for results"
# Creating the directory structure for the Plots so that they dont fill up the main directory and easier to see/mark
if [ ! -d "Plots" ]; then
  mkdir -p Plots
  echo "Created Plots directory"
fi
echo ""
echo ""

#using -O3 for optimization as its a compiler flag for the g++ compiler
# this wont make much of a difference as its not a huge project but still cool to test/know

echo "Running B1"
#Running Part B1 where we calculate the energies and wavefunctions of the hydrogen-like Lithium atom
g++ -O3 hydrogen_like.cpp -o hydrogen -llapack -lblas  && ./hydrogen

echo "Plotting B1"
#Plotting the results (plots stored in Plots directory)
python3 plottingB1.py


echo ""
echo ""


echo "Running B2"
# Running Part B2 where we calculate the energies and wavefunctions of the Lithium atom with a Greens screening potential
g++ -O3 greens.cpp -o greens -llapack -lblas  && ./greens

echo "Plotting B2"
#Plotting the results (plots stored in Plots directory)
python3 plottingB2.py


echo ""
echo ""


echo "Running B3"
# Running Part B3 where we calculate the energies and wavefunctions of the Lithium atom with a 
#Hartree potential using iterative process
g++ -O3 hartree.cpp -o hartree -llapack -lblas  && ./hartree

echo "Plotting B3"
#Plotting the results (plots stored in Plots directory)
python3 plottingB3.py



echo ""
echo ""

echo "Running B4"
# Running Part B4 where we calculate the energies and wavefunctions of the Lithium atom with a 
#Hartree potential using Hartree-Fock iterative method
g++ -O3 hartree-fock.cpp -o hartree-fock -llapack -lblas  && ./hartree-fock

echo "Plotting B4"
#Plotting the results (plots stored in Plots directory)
python3 plottingB4.py



echo ""
echo ""


#Finished!
echo "All calculations and plots are done! Check the Plots directory for the results."
