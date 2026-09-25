#include "LinAlg.h"
#include "Data_structs.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <cmath>

// 1. SpMatBuilder Constructor
SpMatBuilder::SpMatBuilder(int n) : num_rows(n) {}

// 2. Add values into the triplet list
void SpMatBuilder::addVal(int r, int c, double val) {
    triplets.push_back({r, c, val});
}

// 3. Finalize: Sort triplets, sum duplicates, and assemble into CSR format
CSR SpMatBuilder::finalize() {
    CSR mat;
    mat.num_rows = num_rows;
    mat.rowPtr.assign(num_rows + 1, 0);

    if (triplets.empty()) {
        return mat;
    }

    // Sort by row, then by column
    std::sort(triplets.begin(), triplets.end());

    // Merge duplicate contributions at shared nodes
    std::vector<Triplet> merged;
    merged.reserve(triplets.size());

    for (const auto& t : triplets) {
        if (!merged.empty() && merged.back().r == t.r && merged.back().c == t.c) {
            merged.back().val += t.val;
        } else {
            merged.push_back(t);
        }
    }

    // Populate CSR flat arrays, ignoring entries below 1e-15
    mat.vals.reserve(merged.size());
    mat.cols.reserve(merged.size());

    for (const auto& t : merged) {
        if (std::abs(t.val) < 1e-15) {
            continue; // Drop structural/numerical zeros
        }
        mat.vals.push_back(t.val);
        mat.cols.push_back(t.c);
        mat.rowPtr[t.r + 1]++;
    }

    // Prefix sum to compute row offsets
    for (int i = 0; i < num_rows; ++i) {
        mat.rowPtr[i + 1] += mat.rowPtr[i];
    }

    return mat;
}

// 4. Sparse Matrix-Vector Multiplication with dimension validation
void CSR::spmv(const std::vector<double>& x, std::vector<double>& y) const {
    if (x.size() != static_cast<size_t>(num_rows) || y.size() != static_cast<size_t>(num_rows)) {
        throw std::invalid_argument("Vector dimensions do not match CSR matrix.");
    }
    for (int i = 0; i < num_rows; ++i) {
        double sum = 0.0;
        for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j) {
            sum += vals[j] * x[cols[j]];
        }
        y[i] = sum;
    }
}

// 5. Preconditioned Conjugate Gradient Solver with profiling and safety checks
void solve_pcg(const CSR& K, const std::vector<double>& f, std::vector<double>& u, double tol) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int N = K.num_rows;
    if (f.size() != static_cast<size_t>(N)) {
        throw std::invalid_argument("Force vector f dimension does not match matrix rows.");
    }

    // If u was passed in empty or wrong size, allocate and initialize with 0.0
    if (u.size() != static_cast<size_t>(N)) {
        u.assign(N, 0.0);
    }

    // Allocate memory outside the loop (no heap reallocations during solve)
    std::vector<double> r = f; 
    std::vector<double> z(N, 0.0);
    std::vector<double> p(N, 0.0);
    std::vector<double> Ap(N, 0.0);
    std::vector<double> M_inv(N, 1.0);

    // Extract diagonal for Jacobi preconditioner and verify positive-definiteness
    for (int i = 0; i < N; ++i) {
        for (int j = K.rowPtr[i]; j < K.rowPtr[i+1]; ++j) {
            if (K.cols[j] == i) {
                if (K.vals[j] <= 0.0) {
                    throw std::runtime_error("Matrix diagonal is zero/negative. Not Positive-Definite.");
                }
                M_inv[i] = 1.0 / K.vals[j];
                break;
            }
        }
    }

    // Initial residual adjustment
    K.spmv(u, Ap);
    double rz_old = 0.0;
    for (int i = 0; i < N; ++i) {
        r[i] -= Ap[i];
        z[i] = M_inv[i] * r[i];
        p[i] = z[i];
        rz_old += r[i] * z[i];
    }

    int max_iters = N * 2; 
    int iter = 0;
    double max_res = 1.0;

    // Iterative loop
    while (iter < max_iters) {
        K.spmv(p, Ap);
        
        double pAp = 0.0;
        for (int i = 0; i < N; ++i) {
            pAp += p[i] * Ap[i];
        }

        if (pAp <= 0.0) {
            throw std::runtime_error("Matrix indefinite. PCG divergence detected.");
        }

        double alpha = rz_old / pAp;
        double rz_new = 0.0;
        max_res = 0.0;

        for (int i = 0; i < N; ++i) {
            u[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            z[i] = M_inv[i] * r[i];
            rz_new += r[i] * z[i];
            if (std::abs(r[i]) > max_res) {
                max_res = std::abs(r[i]);
            }
        }

        if (max_res < tol) {
            break;
        }

        double beta = rz_new / rz_old;
        for (int i = 0; i < N; ++i) {
            p[i] = z[i] + beta * p[i];
        }
        
        rz_old = rz_new;
        iter++;
    }

    if (iter >= max_iters) {
        throw std::runtime_error("PCG failed to converge within maximum iterations.");
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end_time - start_time;
    
    std::cout << "[PROFILING] PCG Solver Converged in " << iter << " iterations.\n";
    std::cout << "[PROFILING] PCG Execution Time: " << elapsed.count() << " ms.\n";
    std::cout << "[PROFILING] Final Max Residual: " << max_res << "\n";
    std::cout << "[PROFILING] Sparse Matrix Memory Footprint: " << K.get_memory_footprint_mb() << " MB.\n";
}
