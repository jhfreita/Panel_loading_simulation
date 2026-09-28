Orthotropic Panel FEA Solver Based on Mindlin and Timoshenko Plate Theory

Classical Kirchhoff plate theory assumes infinite out-of-plane shear rigidity, a limitation that introduces significant error when analyzing moderately thick structural panels. This implementation utilizes Mindlin-Reissner plate theory to explicitly account for transverse shear deformations and orthotropic behavior.

This project originated from a discrepancy observed throughout my time at Karagozian and case that we did not have enough time to investigate rigorously: verifying the uniform load-deflection relationship for an orthotropic panel with a mixed boundary condition (two fixed edges, one simply supported, and one free edge). An established engineering codebook provided a relationship between aspect ratio (length to width ratio) and stiffness coefficients that was logically inconsistent when compared to other support conditions and implied that adding a fixed edge along the long side of the panel made the panel 5 times less stiff. What began as a manual verification effort scaled into the standalone, object-oriented C++ Finite Element Analysis solver I originally conceptualized during my previous internship but did not have time to build there.

The Engineering Approach

The solver breaks down a continuous wood panel into independent square elements. Because the target material is orthotropic, each element is assigned distinct stiffness capacities in bending, torsion, axial, and shear deformation. The system mathematically assembles these local elements into a global stiffness matrix. 

By applying an external surface load across the system, the engine simultaneously calculates two critical limit states:
1. The maximum physical deflection of the panel.
2. The theoretical maximum load that can be applied before localized internal moments exceed material capacity at the supports. For the scope of this project, shear capacities are neglected, consistent with the scope of my work at Karagozian and Case.

Architecture and Optimization

Global stiffness matrices for highly meshed plates are massive but sparsely populated.

To optimize matrix assembly, I learned about, and implemented, a custom SpMatBuilder function, which uses a Compressed Sparse Row (CSR) format. Instead of allocating a dense 2D array, the builder compresses each row into flat 1D arrays storing non-zero values, column indices, and row pointers indicating which elements of the values vector belong to each row. This eliminates redundant zero-entries, minimizes the memory footprint, and accelerates the Preconditioned Conjugate Gradient (PCG) iterative solver located in LinAlg.cpp.

The broader system architecture was originally designed as one long program but has been updated to break the logic up into physics and formulation:

Physics Formulation (Stiffness.cpp): Computes the 12x12 local stiffness matrix for each element using 2x2 Gauss Quadrature for bending. To prevent shear locking caused by incompatible shape functions (also a new topic I learned about), the engine applies selective reduced 1x1 integration exclusively to the transverse shear terms.

Global Data Structures (Data_structs.h): Manages nodal and elemental configurations, converting complex orthotropic material properties, including independent major and minor bending stiffnesses and rolling shear moduli, into strictly typed C++ structures to prevent memory leaks associated with raw pointers.

Failure Analysis (Failure_check.cpp): Following the PCG displacement calculation, this module sweeps over the elements to map rotations. It calculates localized curvatures at the supports and extracts the internal moment envelopes to pinpoint when failure occurs upon loading.

Build System and Testing

Because FEA relies heavily on floating-point approximations, numerical stability must be continuously verified. The automated test suite (`test_main.cpp`) mathematically validates that the physics engine obeys foundational physical laws by ensuring all generated stiffness matrices maintain Maxwell-Betti symmetry and strict positive-definiteness.

The project utilizes CMake for cross-platform compilation. To build the project from source, execute the following commands in the terminal:

```bash
mkdir build
cd build
cmake ..
make

To run the automated theoretical validation tests, do the following command:
./run_tests

To execute the main finite element solver, do the following command:
./solver
