#pragma once

#include <vector>
#include "solver.hpp"
#include "matrix.hpp"
#include "polynomial.hpp"

template <typename T>
std::vector<Matrix<T>> get_eigenspace(Matrix<T>& mat, T eigen_value) {
    if (!mat.is_square()) throw std::runtime_error("Can't get eigenspace for non-square matrix");
    Matrix<T> mval = Matrix<T>::get_identity(mat.get_rows()) * eigen_value;
    Matrix<T> compute = mat - mval;
    compute.rref();
    return compute.get_non_pivot_columns();
}

template <typename T>
std::vector<T> get_eigenvalues(Matrix<T> mat) {
    if (mat.is_diagonal()) return mat.get_diagonal_elements();

    auto mlambda = Matrix<Polynomial<T>>::get_identity(mat.get_rows()) * Polynomial<T>({{1, 1}}, 'x');
    Matrix<Polynomial<T>> poly_matrix(mat.get_rows(), mat.get_columns());
    std::vector<T> old_flat = mat.get_elements_flat();
    std::vector<Polynomial<T>> new_flat;
    for (const auto x : old_flat) {
        new_flat.push_back(Polynomial(x));
    }
    poly_matrix.init_from_flat(new_flat);

    auto char_matrix = poly_matrix - mlambda;
    auto char_p = char_matrix.get_determinant();

    return solve_polynomial(char_p);
}

template <typename T>
std::pair<Matrix<T>, Matrix<T>> get_diagonalization(Matrix<T>& mat) {
    if (!mat.is_square()) throw std::runtime_error("Can't diagonalize non-square matrix");

    auto eigen_values = get_eigenvalues(mat);
    if (eigen_values.size() != mat.get_columns()) {
        throw std::runtime_error("Matrix is not diagonalizable");
    }

    std::vector<Matrix<T>> eigen_vectors;
    for (const auto eigen_value : eigen_values) {
        auto espace = get_eigenspace(mat, eigen_value);
        for (auto evector : espace) eigen_vectors.push_back(evector);
    }


    Matrix<T> diag(mat.get_rows(), mat.get_columns());
    std::vector<T> diag_flat(eigen_values.size() * eigen_values.size());
    for (int i = 0; i < eigen_values.size(); ++i) {
        for (int j = 0; j < eigen_values.size(); ++j) {
            if (i == j) diag_flat[i * eigen_values.size() + j] = eigen_values[i];
            else diag_flat[i * eigen_values.size() + j] = (T)0;
        }
    }
    diag.init_from_flat(diag_flat);

    Matrix<T> p_matrix(mat.get_rows(), mat.get_columns());
    p_matrix.init_from_columns(eigen_vectors);

    return {diag, p_matrix};
}
