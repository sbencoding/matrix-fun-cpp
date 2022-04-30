#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include <exception>

#define nl '\n'

template <typename T>
class Matrix {

private:
    std::vector<std::vector<T>> data;
    int m;
    int n;

    Matrix<T> get_column(int idx) const {
        if (idx >= n) throw std::invalid_argument("Column index out of bounds");
        Matrix<T> result(n, 1);
        for (int i = 0; i < n; ++i) {
            result.data[i][0] = data[i][idx];
        }
        return result;
    }

    void assert_positive(int num, std::string err_msg) const {
        if (num <= 0) {
            throw std::invalid_argument(err_msg);
        }
    }
public:
    Matrix(const int m, const int n) : m(m), n(n) {
        assert_positive(n, "Matrix dimension has to be positive");
        assert_positive(m, "Matrix dimension has to be positive");
        data = std::vector<std::vector<T>>(m, std::vector<T>(n));
    }

    Matrix(std::vector<std::vector<T>> input) : data(input), m(input.size()), n(input.front().size()) {
        assert_positive(n, "Matrix dimension has to be positive");
        assert_positive(m, "Matrix dimension has to be positive");
    }

    int get_rows() const {
        return m;
    }

    int get_columns() const {
        return n;
    }

    void init_from_flat(const std::vector<T>& flat) {
        if (flat.size() != m * n) throw std::invalid_argument("Number or elements in flat list doesn't match matrix dimension");

        for (int i = 0; i < flat.size(); ++i) {
            data[i / n][i % n] = flat[i];
        }
    }

    std::vector<Matrix<T>> get_column_vectors() const {
        std::vector<Matrix<T>> result;
        for (int i = 0; i < n; ++i) {
            result.push_back(get_column(i));
        }
        return result;
    }

    std::vector<Matrix<T>> get_non_pivot_columns() const {
        // Method assumes matrix is in RREF
        std::vector<Matrix<T>> result;
        int ps_loc = 0;
        for (int i = 0; i < m; ++i) {
            int prev_loc = ps_loc;
            while (ps_loc < n && data[i][ps_loc] != (T)0) ++ps_loc;
            for (int j = prev_loc + 1; j < ps_loc; ++j) {
                result.push_back(get_column(j));
            }
        }
        return result;
    }

    std::vector<T> get_diagonal_elements() const {
        if (!is_square()) throw std::runtime_error("Can't get diagonal elements of non-square matrix");
        std::vector<T> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = data[i][i];
        }
        return result;
    }

    std::vector<T> get_elements_flat() const {
        std::vector<T> flat;
        for (const auto row : data) {
            for (const auto el : row) {
                flat.push_back(el);
            }
        }
        return flat;
    }

    void init_from_columns(const std::vector<Matrix<T>>& cols) {
        if (cols.size() != n) throw std::invalid_argument("Column count mismatch");
        for (int i = 0; i < cols.size(); ++i) {
            for (int j = 0; j < m; ++j) {
                if (cols[i].get_rows() != m || cols[i].get_columns() != 1)
                    throw std::invalid_argument("Dimensions of column vector invalid");
                data[j][i] = cols[i].data[j][0];
            }
        }
    }

    bool is_square() const {
        return m == n;
    }

    bool is_zero() const {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                if (data[i][j] != 0) return false;
            }
        }
        return true;
    }

    bool is_diagonal() const {
        if (!is_square()) return false;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && data[i][j] != 0) return false;
            }
        }
        return true;
    }

    bool is_identity() const {
        if (!is_square()) return false;
        return Matrix<T>::get_identity(n) == *this;
    }

    bool is_lower_triangular() const {
        if (!is_square()) return false;
        for (int i = 0; i < m; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (data[i][j] != 0) return false;
            }
        }
        return true;
    }

    bool is_upper_triangular() const {
        if (!is_square()) return false;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < i; ++j) {
                if (data[i][j] != 0) return false;
            }
        }
        return true;
    }

    Matrix get_inverse() const {
        if (!is_square()) throw std::runtime_error("Unable to get inverse of non-square matrix");

        if (n == 2) {
            Matrix tmp2(data);
            std::swap(tmp2.data[0][0], tmp2.data[1][1]);
            tmp2.data[0][1] *= -1;
            tmp2.data[1][0] *= -1;
            return tmp2 * (1/get_determinant());
        }

        Matrix tmp(m, 2*n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                tmp.data[i][j] = data[i][j];
            }
            for (int j = 0; j < n; ++j) {
                tmp.data[i][j + n] = j == i;
            }
        }
        tmp.rref();

        Matrix result(m, n);
        for (int i = 0; i < m; ++i) {
            for (int j = n; j < 2 * n; ++j) {
                result.data[i][j - n] = tmp.data[i][j];
            }
        }
        return result;
    }

    T get_determinant() const {
        if (!is_square()) throw std::runtime_error("Unable to get determinant of non-square matrix");

        if (n == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }

        T determinant;
        for (int i = 0; i < n; ++i) {
            Matrix tmp(m - 1, n - 1);
            std::vector<T> flat;
            for (int j = 1; j < m; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (k == i) continue;
                    flat.push_back(data[j][k]);
                }
            }
            tmp.init_from_flat(flat);
            determinant += data[0][i] * tmp.get_determinant() * (i & 1 ? -1 : 1);
        }

        return determinant;
    }

    Matrix transpose() const {
        Matrix result(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.data[i][j] = data[j][i];
            }
        }
        return result;
    }

    void rref_sort(std::vector<std::pair<int, std::vector<T>>>& rows) {
        sort(rows.begin(), rows.end(),
            [](const std::pair<int, std::vector<T>>& a, const std::pair<int, std::vector<T>>& b) -> bool {
                if (a.first == b.first) return abs(a.second[a.first]) < abs(b.second[b.first]);
                return a.first < b.first;
            });
    }

    void rref() {
        // Construct leading 0 information
        std::vector<std::pair<int, std::vector<T>>> rows;
        for (auto row : data) {
            int lz = 0;
            while (lz < row.size() && row[lz] == (T)0) ++lz;
            rows.push_back(std::make_pair(lz, row));
        }

        // Get to echelon form
        for (int i = 0; i < m; ++i) {
            // Sort rows according to rref
            rref_sort(rows);
            
            // Index of potential pivot column in current row
            int column = rows[i].first;

            // In case of all zero row we ignore the current row
            if (rows[i].second[column] == 0) continue;

            // Go through the rows below the current one
            for (int j = i + 1; j < m; ++j) {
                // Ignore row if already contains a zero
                if (rows[j].second[column] == 0) continue;

                // Get row multiplier
                T fraction = rows[j].second[column] / rows[i].second[column];
                bool zf = true;
                // Go through each column of the current row
                for (int k = column; k < n; ++k) {
                    // Update column value with fraction
                    rows[j].second[k] -= fraction * rows[i].second[k];
                    
                    // Update leading zero information
                    zf &= (rows[j].second[k] == 0);
                    rows[j].first += zf;
                }
            }
        }

        // Get to row redued echelon form
        for (int i = m - 1; i >= 0; --i) {
            // Get index of potential pivot column
            int column = rows[i].first;
            if (column == n) continue; // ignore zero row
            // Get the value at the pivot location
            T pivot = rows[i].second[column];

            // Divide current row with pivot (so pivot location has a one)
            for (int j = column; j < n; ++j) rows[i].second[j] /= pivot;

            // Go through rows above the current row
            for (int j = i - 1; j >= 0; --j) {
                if (rows[j].second[column] == 0) continue; // ignore row that is already 0 above
                T fraction = rows[j].second[column]; // pivot is already 1 (hence no division by it)
                
                // Update each column of the rows above
                for (int k = column; k < n; ++k) {
                    rows[j].second[k] -= fraction * rows[i].second[k];
                }
            }
        }

        // Change current matrix to the row reduced result
        for (int i = 0; i < m; ++i) {
            data[i] = rows[i].second;
        }
    }

    Matrix pow(int exponent) const {
        Matrix base(data);
        Matrix result(m, n);
        bool rset = false;

        while (exponent) {
            if (exponent & 1) {
                if (!rset) {
                    result.data = base.data;
                    rset = true;
                } else result *= base;
            }
            base *= base;
            exponent >>= 1;
        }
        return result;
    }

    static Matrix get_identity(const int dim) {
        Matrix result(dim, dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                result.data[i][j] = (i == j);
            }
        }
        return result;
    }


    Matrix& operator+=(const Matrix& right) {
        if (right.n != n || right.m != m) throw std::runtime_error("Cannot add matrices with mismatching dimensions");
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                data[i][j] += right.data[i][j];
            }
        }
        return *this;
    }

    friend Matrix operator+(Matrix left, const Matrix& right) {
        left += right;
        return left;
    }

    Matrix& operator-=(const Matrix& right) {
        if (right.n != n || right.m != m) throw std::runtime_error("Cannot subtract matrices with mismatching dimensions");
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                data[i][j] -= right.data[i][j];
            }
        }
        return *this;
    }

    friend Matrix operator-(Matrix left, const Matrix& right) {
        left -= right;
        return left;
    }

    Matrix& operator*=(const T right) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                data[i][j] *= right;
            }
        }
        return *this;
    }

    friend Matrix operator*(Matrix left, const T right) {
        left *= right;
        return left;
    }

    friend Matrix operator*(const T left, Matrix right) {
        right *= left;
        return right;
    }

    Matrix& operator*=(const Matrix& right) {
        if (n != right.m) throw std::runtime_error("#columns on the left matrix doesn't match #rows on the right matrix");
        std::vector<std::vector<T>> result(m, std::vector<T>(n));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < right.n; ++j) {
                for (int k = 0; k < right.m; ++k) {
                    result[i][j] += data[i][k] * right.data[k][j];
                }
            }
        }
        data = result;
        this->n = right.n;
        return *this;
    }

    friend Matrix operator*(Matrix left, const Matrix& right) {
        left *= right;
        return left;
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& obj) {
        os << "[" << nl;
        for (int i = 0; i < obj.m; ++i) {
            os << "  ";
            for (int j = 0; j < obj.n; ++j) {
                os << obj.data[i][j];
                if (j != obj.n - 1) os << ", ";
            }
            os << nl;
        }
        os << "]" << nl;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, Matrix& obj) {
        for (int i = 0; i < obj.m; ++i) {
            for (int j = 0; j < obj.n; ++j) {
                is >> obj.data[i][j];
            }
        }
        return is;
    }

    friend inline bool operator==(const Matrix& left, const Matrix& right) {
        if (left.n != right.n || left.m != right.m) return false;
        return left.data == right.data;
    }
};
