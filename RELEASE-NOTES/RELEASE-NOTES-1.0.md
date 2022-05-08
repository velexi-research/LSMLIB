Release Notes: LSMLIB v1.0
==========================

-------------------------------------------------------------------------------
v1.0.1 (2009/05/03)
===================
* BUG FIXES:
  - Fixed minor typo in INSTALL file. 
  - Fixed function declaration bug in lsm_geometry2d.h.

-------------------------------------------------------------------------------
v1.0.0 (2009/02/15)
===================
* ENHANCEMENTS:
  - Improved support for localization/narrow-band calculations
  - Improved reinitialization via Russo & Smereka subcell fix.
  - Added support for second-order accurate Fast Marching Method 
    calculations:
    * computation of signed distance function
    * solution of Eikonal equation
  - Added support for higher-order delta functions. 
  - Added support for building single precision version of library 
    (via --enable-float configure option).
* MINOR CHANGES:
  - Updated copyright information throughout library.
  - Added many useful helper functions (e.g. curvature calculations, 
    higher-order delta function approximations, etc.).
  - Added support for gfortran compiler (default Fortran compiler) in 
    configuration scripts.
  - Redesigned Fast Marching Method implementation to use generic
    programming ideas.  The field extension, distance function calculation,
    and solution of Eikonal equation are now all "templated" on the
    number of dimensions.
* BUG FIXES:
  - Fixed the following bugs in findLineInTetrahedron():
    * for small tetrahedra, line segment returned occasionally
      extend outside of the specified tetrahedron
    * occationally the intersection point of the line with the
      tetrahedron would be missed if the intersection point lies
      on an edge of the tetrahedron.
  - Fixed several warnings

-------------------------------------------------------------------------------
