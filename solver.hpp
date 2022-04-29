#pragma once

#include <vector>
#include <cmath>
#include "polynomial.hpp"

template<typename T>
std::vector<T> solve_first_order(Polynomial<T> p) {
    const T a = p.get_coefficient(1), b = p.get_coefficient(0);
    return {-b/a};
}

template<typename T>
std::vector<T> solve_second_order(Polynomial<T> p) {
    const T a = p.get_coefficient(2), b = p.get_coefficient(1), c = p.get_coefficient(0);
    const T D = sqrt(b * b -4 * a * c);
    return {
        (-b + D) / (2 * a),
        (-b - D) / (2 * a),
    };
}

template<typename T>
std::vector<T> solve_polynomial(Polynomial<T> polynomial) {
    int deg = polynomial.get_degree();
    switch (deg) {
        case 0: return {polynomial.get_coefficient(0)};
        case 1: return solve_first_order(polynomial);
        case 2: return solve_second_order(polynomial);
    }
    throw std::runtime_error("Can't solver polynomials with degree > 2");
}
