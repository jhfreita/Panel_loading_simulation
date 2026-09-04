#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <array>
#include <algorithm>

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
    
    Material(double EIx, double EIy, double Gx_, double Gy_, double Gxy_, double t_) {
        EIeffx = EIx;
        EIeffy = EIy;
        Gx = Gx_;
        Gy = Gy_;
        Gxy = Gxy_;
        t = t_;
        
        // assumign Poisson effect negligible for CLT; keeping small non zero value 
        nu = 0.02; 
        
        // shear correction factor (rectangular cross section)
        kappa = 5.0 / 6.0; 
        
        // allowable moment capacities (lb-in/in)
        Mx_cap = 10000.0; 
        My_cap = 3000.0; 
    }
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

    void spmv(const std::vector<double>& x, std::vector<double>& y) const {
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

class SpMatBuilder {
public:
   int n;
   std::vector<Triplet> vals_arr;

   SpMatBuilder(int size) {
       n = size;
   }

   void addVal(int r, int c, double val) {
       if (std::abs(val) > 1e-15) {
           vals_arr.push_back({r, c, val});
       }
   }



   CSR finalize() {
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
};

// 12x12 stiffness matrix for Mindlin Q4 element; updated for orthotrnopic CLT
std::array<double, 144> get_element_K(Material mat, double dx, double dy) {
    std::array<double, 144> Ke;
    Ke.fill(0.0);

    double a = dx / 2.0;
    double b = dy / 2.0;
    double t_ = mat.t;

    // constitutive matrix for bending stiffness components
    double Db[3][3] = {0};
    Db[0][0] = mat.EIeffx;
    Db[1][1] = mat.EIeffy;
    
    // geometric mean for Poisson coupling term per standard wood mechanics
    Db[0][1] = mat.nu * std::sqrt(mat.EIeffx * mat.EIeffy); 
    Db[1][0] = Db[0][1];
    Db[2][2] = (mat.Gxy * t_ * t_ * t_) / 12.0; 


    // transverse shear stiffness components using rolling shear moduli
    double Ds[2][2] = {0};
    Ds[0][0] = mat.kappa * mat.Gx * t_;
    Ds[1][1] = mat.kappa * mat.Gy * t_;
    
    
    
    // 2x2 Gauss Quadrature for bending integratio
    double pt = 1.0 / std::sqrt(3.0);
    double g_pts[2] = {-pt, pt};

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
         double xi = g_pts[i];
         double eta = g_pts[j];

         // shape function derivatives w.r.t local coordinates
         double dN_dx[4] = {
             -0.25 * (1.0 - eta) / a,  0.25 * (1.0 - eta) / a,
              0.25 * (1.0 + eta) / a, -0.25 * (1.0 + eta) / a
         };
         
         double dN_dy[4] = {
             -0.25 * (1.0 - xi) / b,  -0.25 * (1.0 + xi) / b,
              0.25 * (1.0 + xi) / b,   0.25 * (1.0 - xi) / b
         };

         // strain-displacement mapping matrix for curvatures
         double B_bend[3][12] = {0};
         for (int k = 0; k < 4; k++) {
             B_bend[0][k * 3 + 1] = dN_dx[k];
             B_bend[1][k * 3 + 2] = dN_dy[k];
             B_bend[2][k * 3 + 1] = dN_dy[k];
             B_bend[2][k * 3 + 2] = dN_dx[k];
         }

         double dV = a * b; 

         for (int m = 0; m < 12; m++) {
             for (int n = 0; n < 12; n++) {
                 double sum = 0.0;
                 for (int p = 0; p < 3; p++) {
                     for (int q = 0; q < 3; q++) {
                         sum += B_bend[p][m] * Db[p][q] * B_bend[q][n];
                     }
                 }
                 Ke[m * 12 + n] += sum * dV;
             }
         }
      }
    }

    // Selective reduced integratio (1x1) for shear terms to preveent locking
    double N0 = 0.25;
    double dN_dx0[4] = {-0.25 / a,  0.25 / a,  0.25 / a, -0.25 / a};
    double dN_dy0[4] = {-0.25 / b, -0.25 / b,  0.25 / b,  0.25 / b};

    double B_shear[2][12] = {0};
    for (int k = 0; k < 4; k++) {
        B_shear[0][k * 3 + 0] = dN_dx0[k];
        B_shear[0][k * 3 + 1] = -N0;
        B_shear[1][k * 3 + 0] = dN_dy0[k];
        B_shear[1][k * 3 + 2] = -N0;
    }

    double dV_shear = 4.0 * a * b; 

    for (int m = 0; m < 12; m++) {
        for (int n = 0; n < 12; n++) {
            double sum = 0.0;
            for (int p = 0; p < 2; p++) {
                sum += B_shear[p][m] * Ds[p][p] * B_shear[p][n];
            }
            Ke[m * 12 + n] += sum * dV_shear;
        }
    }
    
    return Ke;
}



// calc interal momnets at element centroid; compare against orthotropic capacities
void check_capacities(const std::vector<Element>& els, const std::vector<double>& U, Material mat, double dx, double dy, double& mmx, double& mmy) {
    mmx = 0; 
    mmy = 0;
    
    double a = dx / 2.0; 
    double b = dy / 2.0;
    
    double dN_dx[4] = {-0.25/a, 0.25/a, 0.25/a, -0.25/a};
    double dN_dy[4] = {-0.25/b, -0.25/b, 0.25/b, 0.25/b};
    
    double D12 = mat.nu * std::sqrt(mat.EIeffx * mat.EIeffy);
    
    for(size_t i = 0; i < els.size(); i++) {
        double kx = 0; 
        double ky = 0;
        
        for(int j = 0; j < 4; j++) {
            int n_idx = els[i].nodes[j];
            double tx = U[n_idx * 3 + 1];
            double ty = U[n_idx * 3 + 2];
            
            kx += dN_dx[j] * tx; 
            ky += dN_dy[j] * ty; 
        }
        
        // orthotrnpic momnet calculation
        double mx = (mat.EIeffx * kx) + (D12 * ky);
        double my = (D12 * kx) + (mat.EIeffy * ky);
        
        if(std::abs(mx) > mmx) mmx = std::abs(mx);
        if(std::abs(my) > mmy) mmy = std::abs(my);
    }
}


// Preconditioned Conjugate Gradient iterative solver
void solve_pcg(const CSR& A, const std::vector<double>& b, std::vector<double>& x, double tol = 1e-9) {
    int n = b.size();
    x.assign(n, 0.0);

    // Jacobi preconditioner (inverse of diagonal) to accelerate convergence
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

int main() {
    
    // adjust input parameters here; dimensions in inches
    int nx = 21;
    int ny = 21;
    double L = 96.0;
    double W = 96.0;
    
    // initial test load (lbs)
    double base_load = 100.0; 
    
    // CLT panel properties; effective bending stiffnesses, rolling shear moduli, plate thickness
    double input_EIeffx = 115e6;
    double input_EIeffy = 30e6;
    double input_Gx = 50000.0;
    double input_Gy = 15000.0;
    double input_Gxy = 40000.0;
    double input_thickness = 4.125;
    
    Material mat(input_EIeffx, input_EIeffy, input_Gx, input_Gy, input_Gxy, input_thickness);

    double dx = L / (nx - 1);
    double dy = W / (ny - 1);


    // generate nodes across domain sequentially
    std::vector<Node> mesh;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            Node n;
            n.x = i * dx;
            n.y = j * dy;

            // enforce fully clamped boundary conditions along perimeter edges
            if (j == 0 || j == ny - 1 || i == 0 || i == nx - 1) {
                n.fix_w = true; 
                n.fix_tx = true; 
                n.fix_ty = true;
            }
            mesh.push_back(n);
        }
    }

    // build element connectivity matrix; associate four nodes to each quad
    std::vector<Element> els;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            Element el;
            el.nodes[0] = j * nx + i;
            el.nodes[1] = el.nodes[0] + 1;
            el.nodes[2] = (j + 1) * nx + (i + 1);
            el.nodes[3] = (j + 1) * nx + i;
            els.push_back(el);
        }
    }


    int tdofs = mesh.size() * 3;
    SpMatBuilder K_builder(tdofs);
    
    std::vector<double> F(tdofs, 0.0);

    // apply single point load directly to center node
    int mid_node = (ny / 2) * nx + (nx / 2);
    F[mid_node * 3 + 0] = -base_load; 

    // boolean mask to quickly identify constrained degrees of freedom
    std::vector<bool> is_fixed(tdofs, false);
    for (size_t i = 0; i < mesh.size(); i++) {
        if (mesh[i].fix_w)  is_fixed[i * 3 + 0] = true;
        if (mesh[i].fix_tx) is_fixed[i * 3 + 1] = true;
        if (mesh[i].fix_ty) is_fixed[i * 3 + 2] = true;
    }

    std::array<double, 144> Ke = get_element_K(mat, dx, dy);

    // assemble global stiffness matrix
    for (size_t e = 0; e < els.size(); e++) {
      for (int i = 0; i < 4; i++) {
         for (int j = 0; j < 4; j++) {
            for (int di = 0; di < 3; di++) {
               for (int dj = 0; dj < 3; dj++) {
                   int global_i = els[e].nodes[i] * 3 + di;
                   int global_j = els[e].nodes[j] * 3 + dj;

                   // strictly avoid adding values if degree of freedom is fixed; keeps matrix well conditioned
                   if (is_fixed[global_i] || is_fixed[global_j]) continue;

                   int local_i = i * 3 + di;
                   int local_j = j * 3 + dj;
                   K_builder.addVal(global_i, global_j, Ke[local_i * 12 + local_j]);
               }
            }
         }
      }
    }



    // enforce boundary conditios by placing 1.0 on diagonal for all fixed indices
    for (int i = 0; i < tdofs; i++) {
        if (is_fixed[i]) {
            K_builder.addVal(i, i, 1.0);
            F[i] = 0.0;
        }
    }

    CSR K = K_builder.finalize();
    std::vector<double> U_base;
    
    std::cout << "solving with PCG..." << std::endl;
    solve_pcg(K, F, U_base);

    double max_base_d = 0.0;
    for (size_t i = 0; i < mesh.size(); i++) {
        if (std::abs(U_base[i * 3 + 0]) > max_base_d) {
            max_base_d = std::abs(U_base[i * 3 + 0]);
        }
    }
    
    double panel_stiffness = base_load / max_base_d; 
    std::cout << "Panel stiffness: " << panel_stiffness << " lb/in" << std::endl;

    
    double c_ld = 0.0;
    double ld_stp = 500.0; 
    bool failed = false;
    
    std::cout << "\nstepping load..." << std::endl;
    
    
    
    // increment load until internal forces exceed panel material capacity
    while (!failed) {
        c_ld += ld_stp;
        
        // linear analysis; simply scale base deflections by current load ratio
        double scale = c_ld / base_load; 
        
        std::vector<double> U_current(tdofs);
        for(int i = 0; i < tdofs; i++) {
            U_current[i] = U_base[i] * scale;
        }
        
        double mmx;
        double mmy;
        
        check_capacities(els, U_current, mat, dx, dy, mmx, mmy);
        
        
        // Check if maximumn momnet exceeds allowable capacity
        if (mmx >= mat.Mx_cap || mmy >= mat.My_cap) {
            failed = true;
             std::cout << "== PLATE FAILED ==" << std::endl;
             std::cout << "load: " << c_ld / 1000.0 << " kips" << std::endl;
             std::cout << "Momnet X: " << mmx / 1000.0 << " k-in/in" << std::endl;
             std::cout << "Momnet Y: " << mmy / 1000.0 << " k-in/in" << std::endl;
            
            // write final displacement state to file for Paraview visualization
             std::ofstream out("failure_state.vtk");
             out << "# vtk DataFile Version 3.0\n";
             out << "Mindlin Plate\n";
             out << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
             out << "POINTS " << mesh.size() << " double\n";
            
              for (size_t i = 0; i < mesh.size(); i++) {
                 out << mesh[i].x << " " << mesh[i].y << " 0.0\n";
              }
            
             out << "\nCELLS " << els.size() << " " << els.size() * 5 << "\n";
               for (size_t i = 0; i < els.size(); i++) {
                 out << "4 " << els[i].nodes[0] << " " << els[i].nodes[1] << " " << els[i].nodes[2] << " " << els[i].nodes[3] << "\n";
               }
            
              out << "\nCELL_TYPES " << els.size() << "\n";
              for (size_t i = 0; i < els.size(); i++) {
                  out << "9\n";
              }
            
             out << "\nPOINT_DATA " << mesh.size() << "\n";
             out << "SCALARS defelction double 1\nLOOKUP_TABLE default\n";
            
              for (size_t i = 0; i < mesh.size(); i++) {
                 out << U_current[i * 3 + 0] << "\n";
              }
              out.close();
        }
    }
    
    return 0;
}