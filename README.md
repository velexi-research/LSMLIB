Level Set Method Library (LSMLIB)
=================================

-------------------------------------------------------------------------------------------

Table of Contents
-----------------

1. [Overview][#1]

   1.1. [Package Contents][#1.1]

   1.2. [Software Dependencies][#1.2]

   1.3. [License][#1.3]

2. [Installation][#2]

   2.1. [Building the Package][#2.1]

   2.2. [Running Tests][#2.2]

   2.3. [Installing the Package][#2.3]

   2.4. [Uninstalling the Package][#2.4]

3. [Known Issues][#3]

4. [Acknowledgments][#4]

-------------------------------------------------------------------------------------------

## 1. Overview

The Level Set Method Library (LSMLIB) provides support for simulation of implicit surface
and curve dynamics in two- and three-dimensions. It contains an implementation of the
basic level set method algorithms and numerical kernels described in "Level Set Methods
and Dynamics Implicit Surfaces" by S. Osher and R. Fedkiw and "Level Set Methods and Fast
Marching Methods" by J.A. Sethian. It also contains implementations of several advanced
level set method techniques available in the literature.

The library consists of a collection of Fortran subroutines and C/C++ functions.  The main
features of the library include:

* high-computational and high-programmer performance through the use of mixed-language
  programming (e.g. C++/Fortran77, C/Fortran77);

* support for narrow-band/localized computation.

### 1.1. Package Contents

LSMLIB contains a collection of numerical kernels that are commonly used in level set
method calculations. It currently provides numerical kernels for the following:

* computation of spatial derivatives using high-order ENO/WENO schemes

* TVD Runge-Kutta time integration schemes

* imposing boundary conditions

* computation of a variety of geometric quantities (e.g. normal, area/volume, curvature,
  etc.)

* integration over regions defined by implicit functions

* computing the right-hand side of specific level set method partial differential
  equations;

* localization/narrow-band calculations

* fast marching method calculations (e.g. signed distance function, solution to Eikonal
  equation); and

* mathematical functions (e.g. delta functions, norms, etc.).

In addition, LSMLIB provides several utility functions for

* managing of computational grids and data,

* initializing level set functions for simple geometries,

* imposing boundary conditions,

* computing signed distance functions and extension fields via the fast marching method,

* carrying out common data manipulations (e.g. copying data, computing the norm of a grid
  function).

### 1.2. Software Dependencies

* Compilers: C, C++, Fortran

### 1.3. License

See the LICENSE file for copyright and license information.

-------------------------------------------------------------------------------------------

## 2. Installation

### 2.1. Building the Package

* Create a `build` directory, and change into it.

  ```shell
  $ mkdir build
  $ cd build
  ```

* Generate the build files.

  ```shell
  $ cmake ..
  ```

* Build the package.

  ```shell
  $ make
  ```

### 2.2. Running Tests

* From the `build` directory, use `make tests` to build the unit tests.

  ```shell
  $ make tests
  ```

* From the `build` directory, use `ctest` or `make test` to run the unit tests.

  ```shell
  $ ctest
  ```

  For a more detailed test report, use `ctest --verbose`.

  ```shell
  $ ctest --verbose
  ```

### 2.3. Installing the Package

* Set the installation location for the package.

  ```shell
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=/PATH/TO/INSTALL/LOCATION ..
  ```

  __Note__: the `-DCMAKE_INSTALL_PREFIX` option could have been included in
  the command used when generating the build files.

* Use `make install` to install the package.

  ```shell
  $ make install
  ```

* (OPTIONAL) Remove the `build` directory.

### 2.4. Uninstalling the Package

* To uninstall the package, use `make uninstall` from within the `build`
  directory.

  ```shell
  $ cd build
  $ make uninstall
  ```

-------------------------------------------------------------------------------------------

## 3. Known Issues

* Several of the numerical kernels lack unit tests.

* The package only supports one-, two-, and three-dimensional calculations.

* Adapative mesh refinement is not yet available.

* Computation of the maximum stable time step is incorrect for problems involving both
  advection and motion in the normal direction.

* There are a few warnings that show up when using strict compilers (e.g. Intel Compilers).
  These are primarily a result of not explicitly dealing with return values from standard
  C library functions and do not lead to any numerical inaccuracies.

* We do not officially support Windows at this time.

-------------------------------------------------------------------------------------------

## 4. Acknowledgments

The LSMLIB developers gratefully acknowledge support from the following funding agencies:

* Department of Energy under contract numbers DE-FG02-97ER25308 (Computational Science
  Graduate Fellowship) and DE-FC26-06NT43067

* Department of Agriculture under grant #2007-35102-18162

* National Science Foundation

* Air Force Office of Scientific Research

We would also like to thank the following individuals:

* Professor David J. Srolovitz (Yeshiva University) and Professor Steven L. Bryant
  (University of Texas at Austin) - support and encouragement in developing LSMLIB

* Zhaoxuan Wu - original autoconf system for building LSMLIB

* Xiaohai Wan, Thomas Pintelon, and Ning Zhang - documentation for building LSMLIB on
  Windows

* Contributing Users/Bug Reporters: Daniel Thorpe, Markus Gross, Anatoliy Kats, Yi Li,
  Guillaume Walck and his students, Burak Ozkalayci, Zhang Ning, Stefan Sokoll, Moslem
  Kazemi, Danping Zou, Marc Day, and Ruhollah Tavakoli.

* Our humblest apologies if we have accidentally left anybody off of this list. Please let
  us know and we will remedy the situation immediately.

-------------------------------------------------------------------------------------------

[------------------------------------INTERNAL LINKS------------------------------------]: #

[#1]: #1-overview
[#1.1]: #11-package-contents
[#1.2]: #12-software-dependencies
[#1.3]: #13-license

[#2]: #2-installation
[#2.1]: #21-building-the-package
[#2.2]: #22-running-tests
[#2.3]: #23-installing-the-package
[#2.4]: #24-uninstalling-the-package

[#3]: #3-known-issues

[#4]: #4-acknowledgments

[------------------------------------EXTERNAL LINKS------------------------------------]: #
