============================================================================

Level Set Method Library (LSMLIB)
---------------------------------

Kevin T. Chu
Serendipity Research

and

Masa Prodanovic
Institute for Computational Engineering and Sciences, Univ. Texas Austin


============================================================================

Overview
--------

The Level Set Method Library (LSMLIB) provides support for the serial and 
parallel simulation of implicit surface and curve dynamics in two- and 
three-dimensions (support for nodes in one-dimension is also provided). It 
contains an implementation of many of the standard level set method 
algorithms and numerical kernels described in "Level Set Methods and 
Dynamics Implicit Surfaces" by S. Osher and R. Fedkiw and "Level Set Methods 
and Fast Marching Methods" by J.A. Sethian.

The library consists of a collection of Fortran subroutines, C functions, and 
C++ classes. In addition, MATLAB interfaces are provided for some of the basic 
components of the level set method algorithm. The main features of the library 
include:

  * High-computational and high-programmer performance through the use of 
    mixed-language programming (e.g. C++/Fortran77, MATLAB/C/Fortran77);
  * Support for serial and parallel computation;
  * Support for vector level set method calculations;
  * Support for narrow-band; 
  * Support for adaptive mesh refinement (NOT YET IMPLEMENTED);
  * Restart capabilities.

Software Dependencies
---------------------

  * Serial LSMLIB Package: none
  * Parallel LSMLIB Package: SAMRAI (v2.1.0), HDF5 (> v1.6.4), MPI
  * LSMLIB MATLAB Package: MATLAB MEX compiler 

IMPORTANT NOTES:
----------------

  * This library is still under development and has not been thoroughly 
    tested. If you find a bug or have a suggestion, please send me an email 
    at ktchu AT serendipityresearch DOT org, and I will make an attempt to 
    incorporate the appropriate change into the next release of the library. 
    Thank you for your patience during the development process.

  * We are happy to have other researchers use the LSMLIB software in their 
    work (either through collaborations or independently). We merely ask that 
    you keep us aware of your accomplishments and that you reference LSMLIB 
    in presentations and documents that result from your use of the library. 

  * See doc/lsmlib-dox/html/index.html for further documentation on the 
    library.


============================================================================

Contents of Library
-------------------

README.txt ............. this file
INSTALL.txt ............ installation instructions
LICENSE ................ license for LSMLIB
CHANGE_LOG ............. list of changes from version to version
src/.................... directory containing source code
config/................. directory containing configuration files
doc/ ................... directory containing documentation for library
doc/lsmlib-dox/html .... directory containing HTML documentation for library
examples/ .............. directory containing many example programs


============================================================================

Serial LSMLIB Package
---------------------

The Serial LSMLIB package is intended to provide support for serial level set 
method calculations. Currently, it contains functions for

  * managing of computational grids and data,
  * initializing level set functions for simple geometries,
  * imposing boundary conditions,
  * computing signed distance functions and extension fields via the fast 
    marching method, and
  * carrying out common data manipulations (e.g. copying data, computing the 
    max norm of a grid function). 

This package is still actively under development. In the near future, we hope 
to extend the package to cover more of core numerical kernels that support 
level set method calculations:

  * time evolution of the level set function,
  * field extension,
  * reinitialization and orthogonalization (for codimension-two problems), and
  * motion in restricted domains. 


============================================================================

Parallel LSMLIB Package
-----------------------

The Parallel LSMLIB package is built using the Structured Adaptive Mesh 
Refinement Application Infrastructure (SAMRAI) developed and maintained in 
the Center for Applied Scientific Computing (CASC) at Lawrence Livermore 
National (LLNL) The decision to leverage SAMRAI for the Level Set Method 
Library was made because SAMRAI was specifically designed to address several 
key issues:

  * manage parallelism,
  * handle structured, adaptive mesh refinement, and
  * provide check-point and restart capabilities. 

A major design goal of SAMRAI is to provide support for each of these 
features in a scalable manner. Furthermore, the API for SAMRAI is designed 
to shield the application developer from the tedious bookkeeping and memory 
management involved in writing code for parallel, structured adaptive mesh 
refinement applications. The underlying philosophy of SAMRAI is that an 
application developer should be able to migrate to a high-performance, 
parallel program with adaptive mesh capabilities with minimal effort once a 
satisfactory numerical algorithm has been developed for calculations on a 
single, simple rectangular mesh. We have taken advantage of these features 
to implement a high-performance, parallel (eventually adaptive) level set 
method library. 


============================================================================

LSMLIB MATLAB Package
---------------------

The MATLAB LSMLIB package provides a set of MATLAB scripts and MEX-files 
that support level set method calculations. This package currently provides 
functions for

  * time evolution of level set functions,
  * computation of high-order spatial derivatives,
  * total variation diminishing Runge-Kutta time integration, and
  * computation of signed distance functions and extension fields using fast 
    marching methods. 


============================================================================

Acknowledgments
---------------

This software library was supported in part by the Department of Energy 
Computational Science Graduate Fellowship Program of the Office of Science 
and National Nuclear Security Administration in the Department of Energy 
under contract DE-FG02-97ER25308, the National Science Foundation, and the 
Air Force Office for Scientific Research.


============================================================================


Thank you for your interest in LSMLIB!


============================================================================
