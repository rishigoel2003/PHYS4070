#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

#include <stdexcept>

// Forward declarations for LAPACK functions
extern "C" {
    void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, 
               double* w, double* work, int* lwork, int* info);
    void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda,
               double* wr, double* wi, double* vl, int* ldvl, 
               double* vr, int* ldvr, double* work, int* lwork, int* info);
}
using namespace std;
using Complex = complex<double>;
using ComplexVector = vector<Complex>;



class Matrix {
    private:
        std::vector<std::vector<double>> data;
        size_t rows, cols;
    
    public:
        Matrix(size_t r, size_t c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}
        
        Matrix(std::vector<std::vector<double>> values) : data(values), rows(values.size()), cols(values[0].size()) {}
        
        double& operator()(size_t i, size_t j) { return data[i][j]; }
        double operator()(size_t i, size_t j) const { return data[i][j]; }
        
        Matrix tensorProduct(const Matrix& other) const {
            size_t newRows = rows * other.rows;
            size_t newCols = cols * other.cols;
            Matrix result(newRows, newCols);
            
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    for (size_t k = 0; k < other.rows; k++) {
                        for (size_t l = 0; l < other.cols; l++) {
                            result(i * other.rows + k, j * other.cols + l) = data[i][j] * other(k, l);
                        }
                    }
                }
            }
            return result;
        }
        
        void print() const {
            for (const auto& row : data) {
                for (double val : row) {
                    std::cout << val << " ";
                }
                std::cout << "\n";
            }
        }

        // Method to add another matrix (simple element-wise addition)
        Matrix& operator+=(const Matrix& other) {
            for (int i = 0; i < data.size(); ++i) {
                for (int j = 0; j < data[i].size(); ++j) {
                    data[i][j] += other.data[i][j];
                }
            }
            return *this;
        }

        

        // Calculate eigenvalues and eigenvectors using LAPACK
        // Returns a pair: {eigenvalues, eigenvectors}
        // Eigenvectors are stored as columns in the returned matrix
        std::pair<std::vector<double>, Matrix> eigenDecomposition() const {
            // First, ensure the matrix is square
            if (rows != cols) {
                throw std::runtime_error("Eigendecomposition requires a square matrix");
            }

            // LAPACK requires column-major format, so we need to convert our data
            std::vector<double> a(rows * cols);
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    a[j * rows + i] = data[i][j]; // Column-major order
                }
            }

            // Prepare variables for LAPACK's DSYEV routine
            char jobz = 'V'; // Compute both eigenvalues and eigenvectors
            char uplo = 'U'; // Upper triangle of the matrix is stored
            int n = static_cast<int>(rows);
            int lda = n;
            std::vector<double> w(n); // Will hold eigenvalues
            int lwork = 3 * n - 1; // Workspace size
            std::vector<double> work(lwork);
            int info;

            // Call LAPACK's DSYEV routine
            dsyev_(&jobz, &uplo, &n, a.data(), &lda, w.data(), work.data(), &lwork, &info);

            if (info != 0) {
                throw std::runtime_error("LAPACK eigendecomposition failed");
            }

            // Convert back to our Matrix format for eigenvectors
            Matrix eigenvectors(rows, cols);
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    eigenvectors(i, j) = a[j * rows + i]; // Convert from column-major back to row-major
                }
            }

            return {w, eigenvectors};
        }

    };



#endif