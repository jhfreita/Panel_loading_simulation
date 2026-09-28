#include "LinAlg.h"
#include "Data_structs.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>

bool compareTriplets(const Triplet& a, const Triplet& b) {
    if (a.r < b.r) {
        return true;
    }
    else if (a.r > b.r) {
        return false;
    }

    else if (a.c < b.c) {
        return true;
    }
    else if (a.c > b.c) {
        return false;
    }
    else {return false;}
}

SpMatBuilder::SpMatBuilder(int n) {
    num_rows = n;
}

void SpMatBuilder::addVal(int r, int c, double val) {
    Triplet temp_triplet;
    temp_triplet.r = r;
    temp_triplet.c = c;
    temp_triplet.val = val;

    triplets.push_back(temp_triplet);
}

CSR SpMatBuilder::finalize() {
    CSR mat;
    mat.num_rows = num_rows;

    for (int i = 0; i <= num_rows; i = i + 1) {
        mat.rowPtr.push_back(0);
    }

    if (triplets.size() == 0) {
        return mat;
    }

    // Sort using the custom compareTriplets function we defined above
    std::sort(triplets.begin(), triplets.end(), compareTriplets);

    std::vector<Triplet> merged;
    int i = 0;
    
    while (i < (int)triplets.size()) {
        int target_row = triplets[i].r;
        int target_col = triplets[i].c;
        double sum_of_values = 0.0;

        while (i < (int)triplets.size() && triplets[i].r == target_row && triplets[i].c == target_col) {
            sum_of_values = sum_of_values + triplets[i].val;
            i = i + 1; 
        }

        if (std::abs(sum_of_values) > 0.0000001) {
            Triplet combined_triplet;
            combined_triplet.r = target_row;
            combined_triplet.c = target_col;
            combined_triplet.val = sum_of_values;
            merged.push_back(combined_triplet);
        }
    }

    for (int k = 0; k < (int)merged.size(); k++) {
        mat.vals.push_back((merged[k]).val);
        mat.cols.push_back((merged[k]).c);

        int row_index = (merged[k]).r;
        mat.rowPtr[row_index + 1] = mat.rowPtr[row_index + 1] + 1;
    }

    // Prefix sum to convert counts into offsets
    for (int k = 0; k < num_rows; k++) {
        mat.rowPtr[k + 1] = mat.rowPtr[k + 1] + mat.rowPtr[k];
    }

    return mat;
}

void CSR::spmv(const std::vector<double>& x, std::vector<double>& y) const {
    for (int i = 0; i < num_rows; i = i + 1) {
        double row_sum = 0.0;

        int start_index = rowPtr[i];
        int end_index = rowPtr[i + 1];

        for (int j = start_index; j < end_index; j = j + 1) {
            int target_column = cols[j];
            double value = vals[j];

            row_sum = row_sum + (value * x[target_column]);
        }

        y[i] = row_sum;
    }
}

void solve_pcg(const CSR& K, const std::vector<double>& f, std::vector<double>& u, double tol) {
    int N = K.num_rows;

    if (u.size() != (size_t)N) {
        u.clear();
        for (int i = 0; i < N; i = i + 1) {
            u.push_back(0.0);
        }
    }

    std::vector<double> r;      
    std::vector<double> z;      
    std::vector<double> p;      
    std::vector<double> Ap;     
    std::vector<double> M_inv;  

    for (int i = 0; i < N; i = i + 1) {
        r.push_back(f[i]);
        z.push_back(0.0);
        p.push_back(0.0);
        Ap.push_back(0.0);
        M_inv.push_back(1.0);
    }

    for (int i = 0; i < N; i = i + 1) {
        int start_index = K.rowPtr[i];
        int end_index = K.rowPtr[i + 1];

        for (int j = start_index; j < end_index; j = j + 1) {
            if (K.cols[j] == i) {
                double diagonal_val = K.vals[j];

                if (diagonal_val <= 0.0) {
                    throw std::runtime_error("Matrix diagonal is not positive.");
                }

                M_inv[i] = 1.0 / diagonal_val;
                break; 
            }
        }
    }

    K.spmv(u, Ap);

    double rz_old = 0.0;
    for (int i = 0; i < N; i++) {
        r[i] = r[i] - Ap[i];
        z[i] = M_inv[i] * r[i];
        p[i] = z[i];
        rz_old = rz_old + (r[i] * z[i]);
    }

    int max_iterations = N * 2;
    int current_iter = 0;
    double largest_error = 1.0;

    while (current_iter < max_iterations) {
        K.spmv(p, Ap);

        double pAp = 0.0;
        for (int i = 0; i < N; i = i + 1) {
            pAp = pAp + (p[i] * Ap[i]);
        }

        if (pAp <= 0.0) {
            throw std::runtime_error("Matrix is indefinite. PCG cannot continue.");
        }

        double alpha = rz_old / pAp;
        double rz_new = 0.0;
        largest_error = 0.0;

        for (int i = 0; i < N; i = i + 1) {
            u[i] = u[i] + (alpha * p[i]);
            r[i] = r[i] - (alpha * Ap[i]);
            z[i] = M_inv[i] * r[i];

            rz_new = rz_new + (r[i] * z[i]);

            double absolute_r = std::abs(r[i]);
            if (absolute_r > largest_error) {
                largest_error = absolute_r;
            }
        }

        if (largest_error < tol) {
            break;
        }

        double beta = rz_new / rz_old;

        for (int i = 0; i < N; i = i + 1) {
            p[i] = z[i] + (beta * p[i]);
        }

        rz_old = rz_new;
        current_iter = current_iter + 1;
    }

    std::cout << "Solver completed in " << current_iter << " iterations." << std::endl;
    std::cout << "Final max residual error: " << largest_error << std::endl;
}
