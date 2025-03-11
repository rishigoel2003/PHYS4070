#pragma once
#include <cassert>
#include <iostream>
#include <vector>

//==============================================================================

// A relatively basic 'matrix' class. Uses std::vector to store underlying data,
// and provides acces with .at() function, and overload of () operator.
// Of course, this is a very basic example, and has limited features,

class Matrix {

private:
  // use a std::vector to store the data
  std::vector<double> m_data;
  // store the # of rows/cols
  std::size_t m_rows, m_cols;
  // nb: std::size_t is the type of vector.size() - alias for long unsigned int
  // nb: I use 'm_' to signify these are member variables (just convention)

public:
  // contructor: contructs a matrix of dimension: rows X cols
  // by default: matrix will be filled with zero (default from std::vector)
  Matrix(std::size_t rows, std::size_t cols)
      : m_data(rows * cols), m_rows(rows), m_cols(cols) {}

  // Functions are marked 'const' if they do not change class data
  // Return the # of rows
  std::size_t rows() const { return m_rows; }
  // Return the # of columns
  std::size_t cols() const { return m_cols; }

  // defines an 'at' function, access to underlying elements
  // note: this returns a reference to double, so we can modify the data
  double &at(std::size_t i, std::size_t j) {
    assert(i < m_rows && j < m_cols);
    return m_data[i * m_cols + j];
  }
  // a const version of the at() function - returns by copy
  double at(std::size_t i, std::size_t j) const {
    assert(i < m_rows && j < m_cols);
    return m_data[i * m_cols + j];
  }

  // to be fancy, we can also supply a more natural '()' operator to access els
  double &operator()(std::size_t i, std::size_t j) { return at(i, j); }
  double operator()(std::size_t i, std::size_t j) const { return at(i, j); }

  // returns a pointer to first element [i.e., pointer to underlying array]
  // (allows us to use this class as though it were a plain c-style array)
  double *data() {
    return m_data.data();
    // We used the function data() from std::vector here. It is equivilant to:
    // return &m_data[0];
  }
};

//==============================================================================

// 'inline' here just allows us to place the function definitions in the header
// file. If you put these into a .cpp file, it isn't needed.
// (Really, I'm just being lazy)

//----------------------------------------------------------
// overload operators for +=, -=, +, and -
inline Matrix &operator+=(Matrix &a, const Matrix &b) {
  assert(a.rows() == b.rows() && a.cols() == b.cols());
  for (std::size_t i = 0; i < b.rows(); ++i) {
    for (std::size_t j = 0; j < b.cols(); ++j) {
      a(i, j) += b(i, j);
    }
  }
  return a;
}
inline Matrix &operator-=(Matrix &a, const Matrix &b) {
  assert(a.rows() == b.rows() && a.cols() == b.cols());
  for (std::size_t i = 0; i < b.rows(); ++i) {
    for (std::size_t j = 0; j < b.cols(); ++j) {
      a(i, j) -= b(i, j);
    }
  }
  return a;
}
inline Matrix operator+(Matrix a, const Matrix &b) { return a += b; }
inline Matrix operator-(Matrix a, const Matrix &b) { return a -= b; }

//----------------------------------------------------------

// overload operators for scalar multiplication
// allow double*matrix and matrix*double - implement all in terms of 1
inline Matrix &operator*=(Matrix &a, double x) {
  for (std::size_t i = 0; i < a.rows(); ++i) {
    for (std::size_t j = 0; j < a.cols(); ++j) {
      a(i, j) *= x;
    }
  }
  return a;
}
inline Matrix &operator*=(double x, Matrix &a) { return a *= x; }
inline Matrix operator*(Matrix a, double x) { return a *= x; }
inline Matrix operator*(double x, Matrix a) { return a *= x; }
