# FSI-solver
Monolithic 3D solver for Fluid-Structure Interaction written in Matlab and C++ (using Eigen and Spectra). 
The equations are solved in time using a Newmark scheme, during which the nonlinear equations of Large Deformations Elasticity and Fluid Dynamics are solved using an arc-length method.

# Table of Contents
1. [Quick summary](#summary)
2. [Time Stepping](#time)
3. [Nonlinear Solver](#nonlin)
4. [Reduced-Order Modeling](#rom)

<a name="summary"></a>
# Quick summary
This code aims to solve 3D fluid-structure interaction problems consisting of the following coupled equations:
- Large Deformations Elasticity in the solid domain: the materials are assumed to be isotropic and follow a St Venant-Kirchhoff constitutive equation;
- Incompressible Navier-Stokes in the fluid domain: we assume a Newtonian incompressible fluid.

The equations are written using an Arbitrary Lagrangian-Eulerian formulation, and the system is closed by continuity equations for velocity at the interface between fluid and solid, as well as boundary and initial conditions.
A basic method for moving the ALE mesh during the time stepping in order to follow the movement of the interface is implemented. The new mesh is obtained as the solution of an elastic problem defined by displacement boundary conditions assuring that the interface is properly tracked. Thus, distortion is controlled by a pair of elastic coefficients.

The approach being monolithic, a single, large numerical system is solved at each time step. The updating rule at each time increment is given by a Newmark scheme, with adjustable beta and gamma parameters. At each iteration of the Newmark scheme, a nonlinear system of equations has to be solved. 

The solid domain may experience highly nonlinear deformations, and the fluid domain may show a lot of turbulence, which could lead to variants of the Newton-Raphson method failing to converge. Thus we adopt an arc-length method as our nonlinear solver, very robust to overcome strong nonlinearities. The method is more costly than a Newton-Raphson variant but allows here for a single nonlinear solver to solve the nonlinear coupled system directly at each time step, even in the presence of strong nonlinearities.

A method for constructing a reduced-order model using Proper Orthogonal Decomposition will be implemented.

A mesh-generation routine is given in order to build a test-case model with control on the number of elements in the mesh. Block elements are employed here, in a regular 8-node linear form for pressure DOFs and in a 20-node quadratic form for velocity DOFs.

<a name="time"></a>
# Time stepping

<a name="nonlin"></a>
# Nonlinear solver

<a name="rom"></a>
# Reduced-order modeling

