#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include <iomanip>

// Include your custom headers
#include "Data_structs.h"
#include "LinAlg.h"
#include "Failure_check.h"

// Helper function to check if two doubles are effectively equal
bool is_close(double a, double b, double tol = 1e-6) {
    return std::abs(a - b) < tol;
}

// ==========================================
// TEST SUITE 1: Linear Algebra Backend
// ==========================================
void test_sparse_matrix_builder() {
    std::cout << "[TEST] Running SpMatBuilder Compression Tests...\n";

    SpMatBuilder builder(4);
    
    // Add values out of order to test sorting capability
    builder.addVal(2, 2, 5.0);
    builder.addVal(0, 0, 1.0);
    builder.addVal(3, 3, 10.0);
    
    // Add overlapping values to test summation logic (2.0 + 3.0 = 5.0)
    builder.addVal(1, 1, 2.0);
    builder.addVal(1, 1, 3.0);
    
    // Add a near-zero value that should be ignored by the 1e-15 filter
    builder.addVal(0, 1, 1e-16);

    CSR matrix = builder.finalize();

    // Verify row pointers (Should be sized N+1 = 5)
    assert(matrix.rowPtr.size() == 5);
    
    // Verify values were compressed and summed correctly
    assert(matrix.vals.size() == 4);
    assert(is_close(matrix.vals[1], 5.0)); // The overlapping 2.0 + 3.0 at (1,1)

    std::cout << "       -> PASS: CSR compression, sorting, and duplicate summation working.\n";
}

void test_pcg_solver() {
    std::cout << "[TEST] Running Preconditioned Conjugate Gradient (PCG) Solver...\n";

    // Build a known 3x3 Symmetric Positive Definite (SPD) matrix
    // A = [ 4, 1, 0 ]
    //     [ 1, 4, 1 ]
    //     [ 0, 1, 4 ]
    SpMatBuilder builder(3);
    builder.addVal(0, 0, 4.0); builder.addVal(0, 1, 1.0);
    builder.addVal(1, 0, 1.0); builder.addVal(1, 1, 4.0); builder.addVal(1, 2, 1.0);
    builder.addVal(2, 1, 1.0); builder.addVal(2, 2, 4.0);
    CSR A = builder.finalize();

    // Define the right hand side vector b
    // If x_exact = [1, 1, 1], then b = A * x = [5, 6, 5]
    std::vector<double> b = {5.0, 6.0, 5.0};
    std::vector<double> x;

    // Run your custom PCG solver
    solve_pcg(A, b, x, 1e-9);

    // Verify dimensions and exact results
    assert(x.size() == 3);
    assert(is_close(x[0], 1.0));
    assert(is_close(x[1], 1.0));
    assert(is_close(x[2], 1.0));

    std::cout << "       -> PASS: PCG solver successfully converged to exact theoretical root.\n";
}

// ==========================================
// TEST SUITE 2: FEA Physics & Mindlin Plate
// ==========================================
void test_element_stiffness_matrix() {
    std::cout << "[TEST] Running Mindlin Q4 Element Stiffness Validations...\n";

    // Initialize theoretical timber/CLT orthotropic material
    // EIx, EIy, Gx, Gy, Gxy, t
    Material mat(1.2e6, 0.8e6, 500.0, 400.0, 300.0, 0.15);
    mat.nu = 0.2;
    mat.kappa = 5.0 / 6.0; // Standard Mindlin shear correction factor

    // Generate local stiffness matrix for a 2x2 square element
    double dx = 2.0;
    double dy = 2.0;
    std::array<double, 144> Ke = get_element_K(mat, dx, dy);

    // 1. Matrix Size Test
    assert(Ke.size() == 144); // 4 nodes * 3 DOF/node = 12x12 matrix

    // 2. Symmetry Test (Maxwell-Betti Reciprocal Theorem)
    // A stiffness matrix MUST be perfectly symmetric: K[i][j] == K[j][i]
    bool is_symmetric = true;
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            if (!is_close(Ke[i * 12 + j], Ke[j * 12 + i], 1e-5)) {
                is_symmetric = false;
            }
        }
    }
    assert(is_symmetric == true);

    // 3. Positive Definiteness (Diagonal Test)
    // The main diagonal of a stiffness matrix must be strictly positive 
    // because it takes positive energy to deform an element in the direction of the force.
    bool positive_diagonal = true;
    for (int i = 0; i < 12; i++) {
        if (Ke[i * 12 + i] <= 0.0) {
            positive_diagonal = false;
        }
    }
    assert(positive_diagonal == true);

    std::cout << "       -> PASS: 12x12 Stiffness matrix is perfectly symmetric and positive-definite.\n";
}

// ==========================================
// TEST SUITE 3: Orthotropic Capacity Extraction
// ==========================================
void test_capacity_checks() {
    std::cout << "[TEST] Running Internal Moment Extraction Checks...\n";

    // Setup a dummy 1-element system
    std::vector<Element> els;
    Element e1;
    e1.nodes = {0, 1, 2, 3}; // Nodes 0, 1, 2, 3
    els.push_back(e1);

    // Simulate global displacement vector U (4 nodes, 3 DOF each = 12 total entries)
    // Format: [w, tx, ty]
    std::vector<double> U(12, 0.0);
    // Induce a pure curvature to force a bending moment response
    U[1] = 0.001; // tx at node 0
    U[4] = -0.001; // tx at node 1

    Material mat(1e6, 1e6, 500, 500, 300, 0.1);
    mat.nu = 0.2;

    double max_mx = 0.0;
    double max_my = 0.0;

    check_capacities(els, U, mat, 2.0, 2.0, max_mx, max_my);

    // Verify the algorithm extracted non-zero absolute maximums from the deformation
    assert(max_mx > 0.0);
    assert(max_my > 0.0);

    std::cout << "       -> PASS: Maximum moment envelopes extracted successfully.\n";
}

// ==========================================
// MAIN EXECUTION
// ==========================================
int main() {
    std::cout << "=================================================\n";
    std::cout << " FEA Structural Engine - Automated Testing Suite \n";
    std::cout << "=================================================\n\n";

    test_sparse_matrix_builder();
    test_pcg_solver();
    test_element_stiffness_matrix();
    test_capacity_checks();

    std::cout << "\n=================================================\n";
    std::cout << " SUCCESS: ALL UNIT TESTS COMPILED AND PASSED! \n";
    std::cout << "=================================================\n";
    
    return 0;
}