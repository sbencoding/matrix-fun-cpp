#pragma once

#include <map>
#include <cmath>
#include <ostream>

#define nl '\n'

template <typename T>
class Polynomial {
private:
    std::map<int, T> coefficients;
    char variable_letter;
public:
    Polynomial() {
        coefficients[0] = 0;
        variable_letter = 'x';
    }

    Polynomial(const T value) {
        coefficients[0] = value;
        variable_letter = 'x';
    }

    Polynomial(std::map<int, T> coefficients) : coefficients(coefficients) {
        variable_letter = 'x';
    }

    Polynomial(std::map<int, T> coefficients, char variable_letter) : coefficients(coefficients), variable_letter(variable_letter) {}

    bool has_coefficient(int degree) const {
        return coefficients.find(degree) != coefficients.end();
    }

    int get_degree() const {
        if (coefficients.empty()) return 0;
        for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
            if (it->second != (T)0) return it->first;
        }
        return 0;
    }

    T get_coefficient(int degree) const {
        if (!has_coefficient(degree)) return 0;
        return coefficients.at(degree);
    }

    Polynomial& operator+=(const Polynomial& right) {
        for (const auto& kvp : right.coefficients) {
            coefficients[kvp.first] += kvp.second;
        }
        return *this;
    }

    Polynomial& operator+=(const T& right) {
        coefficients[0] += right;
        return *this;
    }

    friend Polynomial operator+(Polynomial left, const Polynomial& right) {
        left += right;
        return left;
    }

    friend Polynomial operator+(T left, Polynomial& right) {
        right += left;
        return right;
    }

    friend Polynomial operator+(Polynomial left, const T& right) {
        left += right;
        return left;
    }

    Polynomial& operator-=(const Polynomial& right) {
        for (const auto& kvp : right.coefficients) {
            coefficients[kvp.first] -= kvp.second;
        }
        return *this;
    }

    Polynomial& operator-=(const T& right) {
        coefficients[0] -= right;
        return *this;
    }

    friend Polynomial operator-(Polynomial left, const Polynomial& right) {
        left -= right;
        return left;
    }

    friend Polynomial operator-(T left, Polynomial& right) {
        right -= left;
        return right;
    }

    friend Polynomial operator-(Polynomial left, const T& right) {
        left -= right;
        return left;
    }

    Polynomial& operator*=(const Polynomial& right) {
        std::map<int, T> tmp;
        for (const auto& outer : coefficients) {
            for (const auto& inner : right.coefficients) {
                tmp[outer.first + inner.first] += outer.second * inner.second;
            }
        }
        coefficients = tmp;
        return *this;
    }

    friend Polynomial operator*(Polynomial left, const Polynomial& right) {
        left *= right;
        return left;
    }

    Polynomial& operator*=(const T& right) {
        for (auto& kvp : coefficients) {
            kvp.second *= right;
        }
        return *this;
    }

    friend Polynomial operator*(T left, Polynomial& right) {
        right *= left;
        return right;
    }

    friend Polynomial operator*(Polynomial left, const T& right) {
        left *= right;
        return left;
    }

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& obj) {
        for (auto it = obj.coefficients.rbegin(); it != obj.coefficients.rend(); ++it) {
            char sign = (it->first >= 0) ? '+' : '-';
            if (it != obj.coefficients.rbegin()) os << ' ' << sign << ' ';
            os << it->second;
            if (it->first != 0) {
                os << " * " << obj.variable_letter << "^" << it->first;
            }
        }
        return os;
    }

    friend inline bool operator==(const Polynomial& left, const Polynomial& right) {
        if (left.coefficients.size() != right.coefficients.size()) return false;
        for (const auto& kvp : left.coefficients) {
            if (!right.has_coefficient(kvp.first)
                    || right.get_coefficient(kvp.first) != kvp.second) return false;
        }
        return true;
    }

    friend inline bool operator!=(const Polynomial& left, const Polynomial& right) {
        return !(left == right);
    }

    Polynomial& operator=(bool value) {
        // Special setter for boolean values to allow behaviour that exists for ints
        // - True is a polynomial with constant factor = 1
        // - False is a polynomail with constant factor = 0
        coefficients[0] = value;
        return *this;
    }

    friend bool operator<(const Polynomial& left, const Polynomial& right) {
        if (left.get_degree() != 0 || right.get_degree() != 0)
            throw std::runtime_error("Can't compare non-constant polynomials");
        return left.get_coefficient(0) < right.get_coefficient(0);
    }

    friend bool operator>(const Polynomial& left, const Polynomial& right) {
        return right < left;
    }

    friend bool operator>=(const Polynomial& left, const Polynomial& right) {
        return !(right > left);
    }

    friend bool operator<=(const Polynomial& left, const Polynomial& right) {
        return !(right < left);
    }
};
