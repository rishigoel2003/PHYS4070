#include "matrix.hpp"
#include "eigen.hpp"
#include <iostream>



int main() {

  //============================================================================
  // Compute eigen{vectors/values} of simple 2x2 matrix:


  // Matrix is a class defined in matrix.hpp
  Matrix sm(2, 2);
  for (int i = 0; i < sm.rows(); ++i) {
    for (int j = 0; j < sm.cols(); ++j) {
      sm(i, j) = 1.0 / (i + j + 1);
    }
  }

  // RealSymmetric: function defined in eigen.hpp
  // MatrixAndVector is a struct defined in eigen.hpp
  // It contians .vector (a std::vector of eigen values)
  // and .matrix (a Matrix of eigenvectors)

  // in c++17, we can do this:
  const auto [EVectors, EValues] =  eigen::RealSymmetric(sm);
  // otherwise, same as this:
  // eigen::MatrixAndVector temp = eigen::RealSymmetric(sm);
  // const auto &EVectors = tmp.matrix;
  // const auto &EValues  = tmp.vector;

  // Print out the eigenvalues and eigenvectors:
  std::cout << "Eigenvalues: ";
  for (auto v : EValues) {
    std::cout << v << ", ";
  }
  std::cout << '\n';
  std::cout << "With corresponding eigenvectors: \n";
  for (int i = 0; i < sm.rows(); ++i) {
    for (int j = 0; j < sm.cols(); ++j) {
      std::cout << EVectors(i, j) << ", ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';

}