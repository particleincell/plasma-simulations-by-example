# About
<img src="https://www.plasmasimulationsbyexample.com/images/plasmasim-book-strip.jpg" alt="Plasma Simulations by Example Cover" width="200" align="right"></img>
This repository contains example codes for Brieda, L., *Plasma Simulations by Example*, CRC Press 2019.
More information about the book can be found on the [companion website](https://www.plasmasimulationsbyexample.com/).

# Plasma Simulations Courses
The material from this book serves as the basis for various plasma simulation classes the author runs on an
annual basis. Current course listing can be found at https://www.particleincell.com/courses/. 

# Organization 
## Chapter 1: Introduction
This chapter begins by introducing the velocity distribution function and differences between fluid and kinetic gas simulation approaches. Next, the basics of particle-based methods are introduced by developing a one-dimensional simulation of an electron oscillating between two grounded electrodes. This section introduces the Leapfrog method, domain discretization, finite difference method (FDM), direct and interactive matrix solvers, and mesh-to-particle interpolation.
## Chapter 2: Grounded Box
The electrostatic particle in cell (ES-PIC) method is introduced by developing a 3D simulation of electrons oscillating in a grounded box filled with background ions. This chapter also serves as a crash course in C++. Topics covered include operator overloading, templates, dynamic memory, and object oriented programming, as well as quiet start, random numbers, and boundary conditions. Simulation results are visualized using Paraview.
## Chapter 3: Flow Around a Sphere
The 3D grounded box ES-PIC code from Chapter 2 is modified to simulate flow of plasma past a stationary sphere (which is analogous to an object moving through through ambient plasma environment). Here we introduce hybrid PIC with Boltzmann relationship for electrons, non-linear Poisson solver, sugarcubing, preconditioned conjugate gradient, steady-state data averaging, calculation of mesh-averaged flow velocity and temperature, and the quasi-neutral approximation.
## Chapter 4: Material Interactions
The example in the prior chapter modeled only ions, which were removed on surface impact to approximate surface neutralization. This chapter modifies the code to explicitly simulate the neutral species. The chapter also reviews velocity sampling from the Maxwellian velocity distribution using uniform and variable macroparticle weight techniques with particle merging. Collisions are introduced next using the Monte Carlo Collisions (MCC) and Direct Simulation Monte Carlo (DSMC) methods. The chapter closes with a discussion of dielectric boundaries and the use of tessellated surface geometries.
## Chapter 5: Symmetry
This chapter introduces approaches for reducing simulation run time by exploiting inherent symmetry or uniformities along a spatial direction. In the case of the flow around the sphere example, symmetry allows one to simulate just a quarter domain. Periodic boundary conditions are also reviewed. Next, a simulation of flow around an infinitely long cylinder is introduced, yielding a planar (XY) 2D simulation. The chapter closes by covering the axisymmetric formulation applicable for simulating cylindrical domains.
## Chapter 6: Unstructured Meshes
This chapter introduces approaches for performing PIC simulations on unstructured meshes. It begins by describing steps needed to generate and load such meshes. Next, a particle tracking code is developed in order to demonstrate how to push particles and interpolate mesh data with unstructured mesh. The finite element method (FEM) is introduced next to produce a Poisson solver. An unstructured version of the flow around a sphere example is then developed.
## Chapter 7: Electromagnetics
This chapter begins by introducing the Boris algorithm for advancing velocities of magnetized particles. Magnetostatic field solver is developed next. This leads to an electromagnetic PIC code utilizing staggered (Yee) grids. The computation of curl, as well as relativistic push are covered.
## Chapter 8: Eulerian Methods
This chapter introduces numerical methods that are based around solving partial differential equations (PDEs) on a stationary grid. The chapter introduces governing equations of magnetohydrodynamics (MHD), with solution methods demonstrated using the model advection-diffusion equation. Stability analysis is covered. The chapter also introduces Vlasov solvers for solving the Boltzmann equation governing the evolution of the velocity distribution function.
## Chapter 9: Parallel Programming
The final chapter introduces parallelization techniques. The chapter covers profiling, multithreading, domain decomposition with MPI, and graphics card processing using CUDA. Common pitfalls, such as race condition, deadlock, and GPU memory transfer overhead are reviewed. A parallel version of the flow around the sphere is developed.

# Building
These examples are meant to be used along with the textbook. Some examples are standalone ``snippets'' used to demonstrate a
particular functionality. Such examples can be compiled and executed using  
```
$ g++ file_name.cpp
$ ./a.out
```  

or, if you prefer specifying the output file name and including optimization  
```
$ g++ -O2 file_name.cpp -o app_name
$ ./app_name
```

Python snippets can be run using  
```
$ python file_name.py
```

Other examples, such as the various iterations of flow around the sphere, consist of multiple source files that need to be compiled
together as in  
```
$ g++ -O2 *.cpp -o sphere
```  
Simulation results are usually stored in a subfolder `results` which needs to be __created manually__ before the code is run for
the first time. The reason for this is that prior to C++17, there was no platform-independent way to create directories
in C++. C++17 was not yet universally supported at the time of book writing. Hence to run the code the first time, use  
```
$ mkdir results
$ ./sphere
```

# Bugs
It is quite likely there are various bugs in the code. Please use the Issues Tracker to identify them so they can be corrected
in a future book revision. 

# Contact
Dr. Lubos Brieda: lubos.brieda@particleincell.com




