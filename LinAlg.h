#pragma once
#include "Data_structs.h"
#include <vector>

class SpMatBuilder {
public:
   int n;
   std::vector<Triplet> vals_arr;

   SpMatBuilder(int size);
   void addVal(int r, int c, double val);
   CSR finalize();
};

// Preconditioned Conjugate Gradient iterative solver
void solve_pcg(const CSR& A, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-9);