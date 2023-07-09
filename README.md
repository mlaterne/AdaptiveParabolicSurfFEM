# AdaptiveParabolicSurfFEM
This repository contains the Code written for a Masters Thesis on Adaptivity for surface FEM for Parabolic surface PDEs based on a posteriori error analysis.
In Particular we solve the PDE: -\Delta_Gamma u + d_t u = f for some close surface \Gamma and some initial condition u_0.

The Code is written in MATLAB and is based on Bartels "Texts in Applied Mathematics 64, Numerical Approximation of Partial Differential Equations"

The Adaptive Algorithm (main.m) is based on the structure presented in Chapter 4.4.
For marking and refining the efficient implmentations for RGB and NVB of Funken are used:
NVB -> "Efficient Implementation of Adaptive P1-FEM in Matlab" S. Funken · D. Praetorius · P. Wissgott
RGB -> "A coarsening algorithm on adaptive red-green-blue refined meshes" Stefan A. Funken 1 · Anja Schmidt 

Further the FEM Solver and many other materials and inspiration was given by B. Kovacs.

The key function to use is the RunScript.m which can be used for testing and generating OUTPUTS for all kinds of surfaces, f and initial conditions.

Some Examples are given in the Examples Folder. The examples showcase some intersting phenomena.
