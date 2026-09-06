#include "Data_structs.h"
#include "LinAlg.h"
#include "Failure_check.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

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

            // enforce fully clamped boundary conditios along perimeter edges
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