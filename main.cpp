#include <iostream>
#include "matrix.hpp"
#include "polynomial.hpp"
#include "eigen.hpp"
#include "approx.hpp"

int main() {
    std::cout << "Driver started" << std::endl;

    /*
     * Basic matrix ops
     */
    Matrix<double> mat(4, 4);
    mat.init_from_flat({1,3,5,9,1,3,1,7,4,3,9,7,5,2,0,9});
    std::cout << "\n===Basic Matrix Ops===\n" << std::endl;
    std::cout << mat << std::endl;
    std::cout << "0.5 multiplier:\n" << mat * 0.5 << std::endl;
    std::cout << "Inverse matrix:\n" << mat.get_inverse() << std::endl;
    std::cout << "square? " << mat.is_square() << std::endl;
    std::cout << "diagonal? " << mat.is_diagonal() << std::endl;
    std::cout << "identity? " << mat.is_identity() << std::endl;
    std::cout << "zero? " << mat.is_zero() << std::endl;
    std::cout << "lower triangular? " << mat.is_lower_triangular() << std::endl;
    std::cout << "upper triangular? " << mat.is_upper_triangular() << std::endl;
    std::cout << "5th power:\n" << mat.pow(5) << std::endl;
    std::cout << "Transpose:\n" << mat.transpose() << std::endl;
    std::cout << "Determinant: " << mat.get_determinant() << std::endl;
    mat.rref();
    std::cout << "RREF:\n" << mat << std::endl;

    /*
     * Eigen values and vectors
     */

    std::cout << "\n===Eigen values and vectors===\n" << std::endl;
    Matrix<double> mat2(2,2);
    mat2.init_from_flat({5,4,1,2});
    std::cout << mat2 << std::endl;

    auto eigen_values = get_eigenvalues(mat2);
    std::cout << "Eigenvalues:\n";
    for (const auto& ev : eigen_values) {
        std::cout << ev << ", ";
    }
    std::cout << std::endl;
    std::cout << "\nDiagonalization:\n";
    auto res = get_diagonalization(mat2);
    std::cout << "D-matrix:\n" << res.first << std::endl;
    std::cout << "P-matrix:\n" << res.second << std::endl;

    /*
     * Approximation and Orthogonalization
     */
    std::cout << "\n===Approximation and Orthogonalization===\n" << std::endl;
    std::cout << "Gram-Schmidt process:\n" << std::endl;
    Matrix<double> g1(3,1), g2(3,1), g3(3, 1);
    g1.init_from_flat({1, -1, 1});
    g2.init_from_flat({1, 0, 1});
    g3.init_from_flat({1, 1, 2});
    std::vector<Matrix<double>> vex = {g1, g2, g3};
    for (int i = 0; i < 3; ++i) {
        std::cout << "Vector " << i + 1 << ":\n";
        std::cout << vex[i] << std::endl;
    }

    auto gs_res = gram_schmidt(vex);
    for (int i = 0; i < 3; ++i) {
        std::cout << "Ortho-Vector " << i + 1 << ":\n";
        std::cout << gs_res[i] << std::endl;
    }


    Matrix<double> v1(3,2), v2(3,1);
    v1.init_from_flat({1,1,1,2,1,4});
    v2.init_from_flat({2,3,3});

    std::cout << "Least squares approximation:\n" << std::endl;
    std::cout << "Coefficient matrix: \n" << v1 << std::endl;
    std::cout << "Solution vector: \n" << v2 << std::endl;

    auto lss = least_squares_approximation(v1,v2);
    std::cout << "Least squares solution is :\n" << lss << std::endl;

    return 0;
}
