#include "LinAlg.h"
#include <cmath>
#include <algorithm>
#include <iostream>

SpMatBuilder::SpMatBuilder(int size) {
    n = size;
}

void SpMatBuilder::addVal(int r, int c, double val) {
    if (std::abs(val) > 1e-15) {
        vals_arr.push_back({r, c, val});
    }
}

CSR SpMatBuilder::finalize() {
    CSR m;
    m.num_rows = n;
    m.rowPtr.assign(n + 1, 0);

    if(vals_arr.empty()) return m;

    std::sort(vals_arr.begin(), vals_arr.end());

    int current_row = -1;
    int last_c = -1;
    double current_sum = 0;

    for(size_t i = 0; i < vals_arr.size(); i++) {
        if(vals_arr[i].r != current_row || vals_arr[i].c != last_c) {
            if(current_row != -1) {
                m.cols.push_back(last_c);
                m.vals.push_back(current_sum);
            }
            
            while(current_row < vals_arr[i].r) {
                current_row++;
                m.rowPtr[current_row] = m.vals.size();
            }
            last_c = vals_arr[i].c;
            current_sum = vals_arr[i].val;
        } else {
            current_sum += vals_arr[i].val; 
        }
    }
    
    m.cols.push_back(last_c);
    m.vals.push_back(current_sum);
    
    while(current_row < n) {
        current_row++;
        m.rowPtr[current_row] = m.vals.size();
    }
    
    return m;
}

void solve_pcg(const CSR& A, const std::vector<double>& b, std::vector<double>& x, double tol) {
    int n = b.size();
    x.assign(n, 0.0);

    std::vector<double> M(n, 1.0); 
    for (int i = 0; i < n; i++) {
       for (int j = A.rowPtr[i]; j < A.rowPtr[i + 1]; j++) {
           if (A.cols[j] == i) {
               if (std::abs(A.vals[j]) > 1e-14) M[i] = 1.0 / A.vals[j];
               break;
           }
       }
    }

    std::vector<double> r = b;
    std::vector<double> z(n);
    for (int i = 0; i < n; i++) z[i] = M[i] * r[i];

    std::vector<double> p = z;
    std::vector<double> Ap(n, 0.0);

    auto dot = [](const std::vector<double>& u, const std::vector<double>& v) {
        double s = 0.0;
        for (size_t i = 0; i < u.size(); i++) s += u[i] * v[i];
        return s;
    };

    double rz_old = dot(r, z);
    double norm_b = std::sqrt(dot(b, b));
    if (norm_b < 1e-14) norm_b = 1.0;

    int max_iters = 5000;
    for (int iter = 0; iter < max_iters; iter++) {
        A.spmv(p, Ap);
        double pAp = dot(p, Ap);

        if (pAp <= 0.0) {
            std::cout << "matrix indefinite error in pcg" << std::endl;
            return;
        }

        double alpha = rz_old / pAp;
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double res_norm = std::sqrt(dot(r, r)) / norm_b;
        if (res_norm < tol) {
            return;
        }

        for (int i = 0; i < n; i++) z[i] = M[i] * r[i];
        
        double rz_new = dot(r, z);
        double beta = rz_new / rz_old;

        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }
        
        rz_old = rz_new;
    }
}