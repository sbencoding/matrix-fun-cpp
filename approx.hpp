#pragma once
#include <vector>
#include "matrix.hpp"

template <typename T>
T vm(const Matrix<T>& v, const Matrix<T>& w) {
    std::vector<T> el_v = v.get_elements_flat(), el_w = w.get_elements_flat();
    T c = (T)0;
    for (int i = 0; i < el_v.size(); ++i) {
        c += el_v[i] * el_w[i];
    }
    return c;
}

template <typename T>
Matrix<T> orthogonal_projection(const Matrix<T>& v, const std::vector<Matrix<T>>& subspace) {
    Matrix<T> result = v;
    for (const auto mat : subspace) {
        T coeff = (vm(v, mat) / vm(mat, mat));
        result -= coeff * mat;
    }
    return result;
}

template<typename T>
std::vector<Matrix<T>> gram_schmidt(const std::vector<Matrix<T>>& vectors) {
    std::vector<Matrix<T>> result;
    if (vectors.empty()) return result;
    result.push_back(vectors.front());
    for (int i = 1; i < vectors.size(); ++i) {
        result.push_back(orthogonal_projection(vectors[i], result));
    }
    return result;
}

template<typename T>
Matrix<T> least_squares_approximation(const Matrix<T>& coeffs, const Matrix<T>& solution) {
    auto lhs = coeffs.transpose() * coeffs;
    auto rhs = coeffs.transpose() * solution;
    Matrix<T> result(lhs.get_rows(), lhs.get_columns() + rhs.get_columns());
    std::vector<T> lhs_flat = lhs.get_elements_flat(), rhs_flat = rhs.get_elements_flat();
    std::vector<T> result_flat(lhs.get_rows() * (lhs.get_columns() + rhs.get_columns()));
    for (int i = 0; i < lhs.get_rows(); ++i) {
        int current_row = i * result.get_columns();

        for (int j = 0; j < lhs.get_columns(); ++j) {
            result_flat[current_row + j] = lhs_flat[current_row + j];
        }
        for (int j = 0; j < rhs.get_columns(); ++j) {
            result_flat[current_row + lhs.get_columns() + j] = rhs_flat[current_row + j];
        }
    }
    result.init_from_flat(result_flat);
    result.rref();
    return result;
}
