# LagrangeMesh-Package

## Description
This is the README file for the Mathematica package LagrangeMesh
version 1.0 (2022) contained in the file LagrangeMesh.wl. LagrangeMesh 
realizes numerically the Lagrange Mesh Method (LMM) in Mathematica.
The basic reference of the method is presented in:

              D. Baye, Phys. Rep. 565, 2015.
              
https://doi.org/10.1016/j.physrep.2014.11.006

In turn, the basic reference of the package can be found in

"Solving the One-Dimensional Time-Independent Schrödinger Equation with 
      High Accuracy: The LagrangeMesh Mathematica Package"
                    by J.C. del Valle, 2022
                     jcdvaller@gmail.com

              (https://arxiv.org/abs/2208.14340)



## NNotation

   * E: eigenvalues/energy
   * m: mass of the particle
   * x: real variable
* y(x): eigenfunction/wavefunction in variable x
* V(x): confining potential




LagrangeMesh 1.0 is a Mathematica package devoted to solving numerically 
the one-dimensional time-independent Schrödinger equation with Dirichlet
boundary conditions:

     -(1/2m)*y''(x)+V(x)*y(x)=E*y(x) , y(0)=y(b)=0 ,      a < x < b ,

for an arbitrary real interval (a,b). 

===

(*INSTALLATION*)

0. Download the file LagrangeMesh.wl

A. *Full Installation*. Open Mathematica, go to File -> Install.
   A window will prompt you for more information. Fill in as follows:

  i. Type of Item to Install -> Package
 ii. Source -> Select file LagrangeMesh.wl 
iii. Install Name -> LagrangeMesh
 iv. Click OK.

Once the installation was successful, the file LagrangeMesh.wl can be
deleted. To load the package in a given NoteBook, type Needs["LagrangeMesh`"] 
in the command line and evaluate the cell. Each time the kernel is restarted, 
the package has to be loaded as indicated.


B. *Loading the Package*. This is a simple alternative that requires no
   full installation. Loading the package is achieved by evaluating
   <<"path/LagrangeMesh.wl" in the command line of the notebook
   we want to work with. The explicit form of path is found evaluating
   the command NotebookDirectory[]. This option only works if the package
   and the notebook (in which we will perform calculations) are contained
   in the same directory. Each time the kernel is restarted, package has to be
   loaded as indicated.

===

(*COMMANDS*)

Once installed/loaded, the package provides the user with five commands:

Name                                    Goal

BuildMesh              --   Construction of mesh points and weights            
AvailableMeshQ         --   Check available meshes and weights
LagMeshEigenvalues     --   Computes eigenvalues
LagMeshEigenfunctions  --   Computes eigenfunctions
LagMeshEigensystem     --   Computes eigenvalues and eigenfunctions


===

(*USAGE*)

Based on particular examples, we describe the usage of each command. Input and Output
are shown in Mathematica's style. We assume that all evaluations are performed from the
notebook called MyWorkNotebook.nb that is contained in a folder called MyDirectory. By
default the mass is set m = 1.


(1) BuildMesh. It constructs the requested mesh points and weights for a given type 
               of classical orthogonal polynomial. (Legendre,Laguerre,Hermite).
               Basic example: 

 In[1]:= BuildMesh["Hermite",15,WorkingPrecision->20,Weights->True]
Out[1]:= Hermite_15_WP_20.dat
         Hermite_15_WP_20_Weights.dat 
                
               As Output, shown in Out[1], the program prints on screen the name of
               two files that were generated and stored. Meshes and weights are 
               automatically stored according to the following tree diagram:

 MyDirectory
 |-- MyWorkNotebook.nb
 `-- Meshes
     |-- Hermite
     |   |-- MeshPoints
     |   |   |-- Hermite_15_WP_20.dat
     |   `-- Weights
     |       `-- Hermite_15_WP_20_Weights.dat
     |-- Laguerre 
     `-- Legendre

All directories inside MyDirectory will be automatically created on the first use
of the command BuildMesh.       


(2) AvailableMeshQ. By specifying Type, Dimension, and WorkingPrecision, it delivers True on 
                    screen if the mesh was previously constructed. Otherwise, the Output
                    is False. In addition, it can deliver an ordered table that shows
                    all stored meshes. A basic example:

 In[2]:= AvailableMeshQ["Hermite",WorkingPrecision->20,Dimension->15]
Out[2]:= True
                                 



(3) LagMeshEigenvalues. It computes the eigenvalues for a given potential defined on
                        an interval. The number of mesh points must be specified as 
                        well as the requested lowest eigenvalues. If WorkingPrecision
                        is not specified, calculations are performed with MachinePrecision.
                        A basic example is presented below: the first three eigenvalues
                        of the harmonic oscillator defined on the whole real line
                        calculated with 15 mesh points and WorkinPrecision->20.

 In[3]:= LagMeshEigenvalues[1/2*x^2,{x,-Infinity,Infinity},3,15,WorkingPrecision->20]
Out[3]:= {0.500000000000000000,1.50000000000000000,2.50000000000000000}
          
 

(4) LagMeshEigenfunctions. It computes the eigenfunctions for a given potential defined on
                           an interval. The number of mesh points must be specified as 
                           well as the requested lowest eigenvalues. If WorkingPrecision
                           is not specified, calculations are performed with MachinePrecision.
                           A basic example is presented below: the first eigenfunction 
                           of the harmonic oscillator defined on the whole real line
                           calculated with 15 mesh points and WorkinPrecision->20.

 In[4]:= LagMeshEigenfunctions[1/2*x^2,{x,-Infinity,Infinity},1,15,WorkingPrecision->20]
Out[4]:= -0.75112554446494248283 E^(-(x^2/2))

(5) LagMeshEigensystem. It computes the eigenvalues and eigenfunctions simultaneously 
                        for a given potential defined on an interval. The number of mesh
                        points must be specified as well as the requested lowest eigenvalues.
                        If WorkingPrecision is not specified, calculations are performed with
                        MachinePrecision. A basic example is presented below: the first
                        eigenvalue and eigenfunction of the harmonic oscillator defined
                        on the whole real line calculated with 15 mesh points and
                        WorkinPrecision->20.

 In[5]:= LagMeshEigensystem[1/2*x^2,{x,-Infinity,Infinity},1,15,WorkingPrecision->20]
Out[5]:= {0.500000000000000000, -0.75112554446494248283 E^(-(x^2/2))}

===

LagrangeMesh 1.0 (2022):

(C) J.C. del Valle, Version 1.0 (September,2022)  
    (Written in Mathematica 13 and tested in version 12 and 13)

Please report any malfunction of the package to jcdvaller@gmail.com
