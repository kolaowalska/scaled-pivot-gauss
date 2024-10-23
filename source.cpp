#include "vectalg.h"
#include <cmath>
#include <numeric>

Vector elimination(const Matrix & A, const Vector & b, Vector scale) {
    size_t n = A.size();
    Matrix a(A);
    Vector b_copy(b), x(n), p(n);

    std::iota(p.begin(), p.end(), 0);

    for (size_t k = 0; k < n - 1; ++k) {
        size_t maxRow = k;
        double maxScaled = std::fabs(a(k, k) / scale[k]);

        for (size_t i = k + 1; i < n; ++i) {
            double scaled = std::fabs(a(i, k) / scale[i]);
            if (scaled > maxScaled) {
                maxScaled = scaled;
                maxRow = i;
            }
        }

        std::swap(p[k], p[maxRow]);
        std::swap(scale[k], scale[maxRow]);

        for (size_t i = k + 1; i < n; ++i) {
            double z = a(p[i], k) / a(p[k], k);
            a(p[i], k) = z;
            for (size_t j = k + 1; j < n; ++j) {
                a(p[i], j) -= z * a(p[k], j);
            }
        }
    }

    for (size_t k = 0; k < n - 1; ++k) {
        for (size_t i = k + 1; i < n; ++i) {
            b_copy[p[i]] -= a(p[i], k) * b_copy[p[k]];
        }
    }

    x[n - 1] = b_copy[p[n - 1]] / a(p[n - 1], (n - 1));
    for (int i = n - 2; i >= 0; --i) {
        double sum = 0.0;
        for (size_t j = i + 1; j < n; ++j) {
            sum += a(p[i], j) * x[j];
        }
        x[i] = (b_copy[p[i]] - sum) / a(p[i], i);
    }
    return x;
}

Vector solveEquations(const Matrix & A, const Vector & b, double eps){
    size_t n = A.size();
    Vector scale(n);

    for (int i = 0; i < A.size(); ++i){
        Vector aux(A.size());
        for (int j = 0; j < A.size(); ++j){
            aux[j] = A(i, j);
        }
        scale[i] = aux.max_norm();
    }

    Vector residual = residual_vector(A, b, scale);
    while (true) {
        Vector margin = elimination(A, residual_vector(A, b, residual), scale);
        for (size_t i = 0; i < residual.size(); ++i) {
            residual[i] += margin[i];
        }
        if(residual_vector(A, b, residual).max_norm() < eps) break;
    }
    return residual;
}
