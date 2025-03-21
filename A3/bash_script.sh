set -e  # Exit on error



if [ ! -d "Plots" ]; then
  mkdir -p Plots
  echo "Created Plots directory"
fi


# g++ hydrogen_like.cpp -o hydrogen -llapack -lblas  && ./hydrogen


# python3 plottingB1.py


# g++ greens.cpp -o greens -llapack -lblas  && ./greens


# python3 plottingB2.py



g++ hartree.cpp -o hartree -llapack -lblas  && ./hartree

python3 plottingB3.py