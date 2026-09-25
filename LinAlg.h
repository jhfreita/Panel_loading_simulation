#pragma once
#include "Data_structs.h"
#include <vector>


void solve_pcg(const CSR& K, const std::vector<double>& f, std::vector<double>& u, double tol = 1e-6);
