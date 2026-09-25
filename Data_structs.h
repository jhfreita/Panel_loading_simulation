#pragma once
#include <vector>
#include <array>
#include <stdexcept>

// Assumptions: Symmetric cross section; perfectly orthotropic 
// Neglect P-delta effects; Fourier loading based on Timoshenko
struct Material {
    double EIeffx;
    double EIeffy;
    double Gx;
    double Gy;
    double Gxy;
    
    double nu;
    double t;
    double kappa; 
    
    double Mx_cap; 
    double My_cap;
    
    // Rigorous Constructor: Validates thermodynamic and physical reality
    Material(double EIx, double EIy, double Gx_, double Gy_, double Gxy_, double t_) {
        if (t_ <= 0.0) {
            throw std::invalid_argument("Material thickness must be strictly positive.");
        }
        if (EIx <= 0.0 || EIy <= 0.0) {
            throw std::invalid_argument("Flexural rigidities must be strictly positive.");
        }
        if (Gx_ <= 0.0 || Gy_ <= 0.0 || Gxy_ <= 0.0) {
            throw std::invalid_argument("Shear moduli must be strictly positive.");
        }
        
        EIeffx = EIx;
        EIeffy = EIy;
        Gx = Gx_;
        Gy = Gy_;
        Gxy = Gxy_;
        t = t_;
        
        nu = 0.02; 
        kappa = 5.0 / 6.0;
        Mx_cap = 10000.0; 
        My_cap = 3000.0;
    }
};

struct Node {
    double x;
    double y;
    bool fix_w;
    bool fix_tx;
    bool fix_ty;
    
    Node() : x(0.0), y(0.0), fix_w(false), fix_tx(false), fix_ty(false) {}
};

struct Element {
    std::array<int, 4> nodes;
};

struct CSR {
    int num_rows;
    std::vector<double> vals;
    std::vector<int> cols;
    std::vector<int> rowPtr;

    void spmv(const std::vector<double>& x, std::vector<double>& y) const;
    
    double get_memory_footprint_mb() const {
        size_t bytes = 0;
        bytes += vals.capacity() * sizeof(double);
        bytes += cols.capacity() * sizeof(int);
        bytes += rowPtr.capacity() * sizeof(int);
        return static_cast<double>(bytes) / (1024.0 * 1024.0);
    }
};

struct Triplet {
    int r;
    int c;
    double val;
    
    bool operator<(const Triplet& other) const {
        if (r == other.r) return c < other.c;
        return r < other.r;
    }
};

// The missing builder required for LinAlg.cpp
struct SpMatBuilder {
    int num_rows;
    std::vector<Triplet> triplets;
    
    SpMatBuilder(int n);
    void addVal(int r, int c, double val);
    CSR finalize();
};
