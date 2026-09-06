#pragma once
#include <vector>
#include <array>

// Assumptions: Symmetric cross section; perfectly orthotropic (major and minor stiffnesses 
// aligned perfectly with the x ad y axes)
// Neglect P-delta effects; Fourier loading based on Timoshenko "Theory of Plates and Shells" (1959)
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
    
    Material(double EIx, double EIy, double Gx_, double Gy_, double Gxy_, double t_);
};

// Global coordinate structure; holds physical coordinates adn boundary condition flags
struct Node {
    double x;
    double y;
    bool fix_w;
    bool fix_tx;
    bool fix_ty;
    
    Node() : x(0), y(0), fix_w(false), fix_tx(false), fix_ty(false) {}
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