#include "Failure_check.h"
#include <cmath>

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