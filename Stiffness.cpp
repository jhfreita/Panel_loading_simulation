#include "Data_structs.h"

Material::Material(double EIx, double EIy, double Gx_, double Gy_, double Gxy_, double t_) {
    EIeffx = EIx;
    EIeffy = EIy;
    Gx = Gx_;
    Gy = Gy_;
    Gxy = Gxy_;
    t = t_;
    
    // assumign Poisson effect negligible for CLT; keeping small non zero value 
    nu = 0.02; 
    
    // Mindlin shear correction factor (rectangular cross section)
    kappa = 5.0 / 6.0; 
    
    // allowable moment capacities (lb-in/in)
    Mx_cap = 10000.0; 
    My_cap = 3000.0; 
}

void CSR::spmv(const std::vector<double>& x, std::vector<double>& y) const {
    for (int i = 0; i < num_rows; i++) {
        double dot = 0.0;
        int start = rowPtr[i];
        int end = rowPtr[i + 1];
        for (int j = start; j < end; j++) {
            dot += vals[j] * x[cols[j]];
        }
        y[i] = dot;
    }
}