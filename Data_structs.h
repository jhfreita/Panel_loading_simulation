#pragma once

#include <vector>
#include <array>
#include <cstddef>

struct Material {
    double EIeffx;
    double EIeffy;
    double Gx;
    double Gy;
    double Gxy;
    double thickness;
    double t;
    double Mx_cap;
    double My_cap;
    double nu;     // Poisson's ratio
    double kappa;  // Mindlin shear factor

    Material() {
        EIeffx = 0.0;
        EIeffy = 0.0;
        Gx = 0.0;
        Gy = 0.0;
        Gxy = 0.0;
        thickness = 0.0;
        t = 0.0;
        Mx_cap = 0.0;
        My_cap = 0.0;
        nu = 0.2;
        kappa = 5.0 / 6.0;
    }

    Material(double ex, double ey, double gx, double gy, double gxy, double thick) {
        EIeffx = ex;
        EIeffy = ey;
        Gx = gx;
        Gy = gy;
        Gxy = gxy;
        thickness = thick;
        t = thick;
        Mx_cap = 0.0;
        My_cap = 0.0;
        nu = 0.2;
        kappa = 5.0 / 6.0;
    }
};


struct Node {
    int id;
    double x;
    double y;

    // Boundary condition flags (Degrees of Freedom)
    bool fix_w;   // Out-of-plane displacement
    bool fix_tx;  // Rotation about x
    bool fix_ty;  // Rotation about y

    Node() {
        id = 0;
        x = 0.0;
        y = 0.0;
        fix_w = false;
        fix_tx = false;
        fix_ty = false;
    }
};


struct Element {
    int id;
    std::array<int, 4> nodes;

    Element() {
        id = 0;
        nodes[0] = 0;
        nodes[1] = 0;
        nodes[2] = 0;
        nodes[3] = 0;
    }
};


struct Triplet {
    int r;
    int c;
    double val;

    Triplet() {
        r = 0;
        c = 0;
        val = 0.0;
    }

    Triplet(int row_index, int col_index, double value_in) {
        r = row_index;
        c = col_index;
        val = value_in;
    }
};

bool compareTriplets(const Triplet& a, const Triplet& b);

struct CSR {
    int num_rows;
    std::vector<int> rowPtr;
    std::vector<int> cols;
    std::vector<double> vals;

    CSR() {
        num_rows = 0;
    }

    void spmv(const std::vector<double>& x, std::vector<double>& y) const;

    double get_memory_footprint_mb() const {
        double total_bytes = 0.0;

        double row_ptr_bytes = (double)(rowPtr.size() * sizeof(int));
        total_bytes = total_bytes + row_ptr_bytes;

        double cols_bytes = (double)(cols.size() * sizeof(int));
        total_bytes = total_bytes + cols_bytes;

        double vals_bytes = (double)(vals.size() * sizeof(double));
        total_bytes = total_bytes + vals_bytes;

        double megabytes = total_bytes / (1024.0 * 1024.0);
        return megabytes;
    }
};

struct SpMatBuilder {
    int num_rows;
    std::vector<Triplet> triplets;

    SpMatBuilder(int n);
    void addVal(int r, int c, double val);
    CSR finalize();
};
