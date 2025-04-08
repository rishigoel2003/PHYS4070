#pragma once
#include "matrix.hpp"
#include <vector>
#include <cassert>
#include <iostream>

// Very simple struct; stores a matrix and a vector - used to hold the matrix of
// eigenvectors, and a vector of eigenvalues
struct MatrixAndVector {
  MatrixAndVector(std::size_t dimension)
      : matrix(dimension, dimension), vector(dimension) {}
  Matrix matrix;
  std::vector<double> vector;
};




extern "C"
    // dsyev_ is a symbol in the LAPACK library files
    // Documentation: http://www.netlib.org/lapack/explore-html/index.html
    int
    dsyev_(char *jobz, char *uplo, int *dimension, double *in_out_matrix,
           int *dimension2, double *out_e_values, double *workspace_array,
           int *workspace_size, int *info);
// nb: The variable names in the declaration here are not required, but I find
// them helpful
// This can go either in the .cpp or the .hpp
extern "C"
    int 
    dsygv_(
      int *ITYPE,   // =1 for problems of type Av=eBv
      char *JOBZ,   // ='V' means calculate eigenvectors
      char *UPLO,   // 'U': upper triangle of matrix is stored, 'L': lower
      int *N,       // dimension of matrix A
      double *A,    // c-style array for matrix A (ptr to array, pointer to a[0])
                    // On output, A contains matrix of eigenvectors
      int *LDA,     // For us, LDA=N
      double *B,    // c-style array for matrix B [Av=eBv]
      int *LDB,     // For us, LDB =N
      double *W,    // Array of dimension N - will hold eigenvalues
      double *WORK, // 'workspace': array of dimension LWORK
      int *LWORK,   // dimension of workspace: ~ 6*N works well
      int *INFO     // error code: 0=worked.
   );


//------------------------------------------------------------------------------
// Finds eigen values and eigen vectors of real, symmetric matrix
// Returns a 'MatrixAndVector' - defined in eigen.hpp
// note: We take the matrix by copy. This is because dsyev_ destroys the
// matrix. If we no longer need 'matrix' can call like this:
// RealSymmetric(std::move(matrix));
// This will 'move' the matrix, instead of copying it, and will destroy
// the original matrix in the process.
MatrixAndVector RealSymmetric(Matrix matrix) {
  assert(matrix.rows() == matrix.cols() &&
         "RealSymmetric only works for square matrix");

  MatrixAndVector result(matrix.rows());

  // Data required by dsyev. See documentation for description of each option:
  // http://www.netlib.org/lapack/explore-html/index.html
  char jobz{'V'};
  char uplo{'U'};
  int dimension = static_cast<int>(matrix.rows()); // LAPACK expects an int
  int lwork = 6 * dimension;
  std::vector<double> work(static_cast<std::size_t>(lwork)); //making a vector of size lwork
  int info;
  // calculate eigenvalues using the DSYEV lapack subroutine
  // on INPUT, 'matrix' is the matrix to invert. After dsyev_ exits, 'matrix'
  // will hold the eigenvectors, and 'result.vector' will hold the eigenvalues
  dsyev_(&jobz, &uplo, &dimension, matrix.data(), &dimension,
         result.vector.data(), work.data(), &lwork, &info);

  // move the data from original matrix into matrix that lives in
  // 'MatrixAndVector' - which we return
  result.matrix = std::move(matrix);

  // check for errors
  if (info != 0) {
    std::cout << "D'oh! DSYEV returned error code: " << info << '\n';
  }

  return result; // nb: named return-value optimisation (NRVO) is used here
}













MatrixAndVector GeneralisedEigenvalue(Matrix A, Matrix B) {
  assert(A.rows() == A.cols() &&
         "Generalised Eigenvalue only works for square matrix");
  assert(B.rows() == B.cols() &&
         "Generalised Eigenvalue only works for square matrix");
  assert(A.rows() == B.cols() &&
         "Generalised Eigenvalue only works for square matrix");

  MatrixAndVector result(A.rows());

  // Data required by dsygv. See documentation for description of each option:
  // http://www.netlib.org/lapack/explore-html/index.html
  char jobz{'V'};
  char uplo{'U'};
  int dimension = static_cast<int>(A.rows()); // LAPACK expects an int
  int lwork = 6 * dimension;
  std::vector<double> work(static_cast<std::size_t>(lwork)); //making a vector of size lwork
  int info;
  int type = 1;

  // calculate eigenvalues using the DSYGV lapack subroutine
  // on INPUT, 'matrix' is the matrix to invert. After dsygv_ exits, 'matrix'
  // will hold the eigenvectors, and 'result.vector' will hold the eigenvalues
  dsygv_(&type,&jobz, &uplo, &dimension, A.data(), &dimension, B.data(), &dimension,
         result.vector.data(), work.data(), &lwork, &info);

  // move the data from original matrix into matrix that lives in
  // 'MatrixAndVector' - which we return
  result.matrix = std::move(A);

  // check for errors
  if (info != 0) {
    std::cout << "D'oh! DSYGV returned error code: " << info << '\n';
  }


  // for (int j = 0; j < dimension; ++j) {
  //   double norm = 0.0;
  //   for (int i = 0; i < dimension; ++i) {
  //     norm += result.matrix(i, j) * result.matrix(i, j);
  //   }
  //   norm = std::sqrt(norm);
    
  //   if (norm > 1e-10) {  // Avoid division by very small numbers
  //     for (int i = 0; i < dimension; ++i) {
  //       result.matrix(i, j) /= norm;
  //     }
  //   }
  // }
  


  return result; // nb: named return-value optimisation (NRVO) is used here
}



