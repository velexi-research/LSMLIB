Release Notes: LSMLIB v0.7
==========================

-------------------------------------------------------------------------------
v0.7.2 (2006/03/11)
===================
* ADDED FEATURES:
  - The library now uses the autoconf system for building the library.
  - Modified LSM3D_findLineInTetrahedron() so that it returns the endpoints
    of the line segment {phi=0,psi=0} oriented in the same direction as
    grad(phi) cross grad(psi)
* REMOVED FEATURES:
  - Removed ChomboVisDataWriter class since it is easier to use the VisIt
    data writer capabilities provided by SAMRAI.
* BUG FIXES:
  - Modified computation of unit normals so that they actually have norm 1 
    to machine precision (as opposed to having O(dx) error).
  - Modified computation of stable time step to avoid division by zero if
    the velocity is identically zero throughout the computational domain.
  - GCC 4.0 compatibility
    * Added explicit qualifier for the SAMRAI tbox::MPI struct.
    * Fixed syntax bug in ReinitializationAlgorithm, FieldExtensionAlgorithm
      OrthogonalizationAlgorithm

-------------------------------------------------------------------------------
v0.7.1 (2006/02/10)
===================
* ADDED FEATURES:
  - Added reinitializeLevelSetFunctions() and orthogonalizeLevelSetFunctions()
    to the LevelSetMethodAlgorithm class so that users can manually
    reinitialize and orthogonalize the level set functions.
  - Changed default stopping condition for ReinitializationAlgorithm,
    FieldExtensionAlgorithm, and OrthogonalizationAlgorithm to use
    the stop distance criterion with the stop distance set to be the
    length of the largest dimension of the computational domain.
  - Improved verbose mode status output by ReinitializationAlgorithm,
    FieldExtensionAlgorithm, and OrthogonalizationAlgorithm. 
  - Minor improvements to documentation
* API CHANGES:
  - Renamed LevelSetMethodToolbox<>::copyDataToScratchSpace() to 
    copySAMRAIData().
  - Reordered arguments to constructors of ReinitializationAlgorithm,
    FieldExtensionAlgorithm, and OrthogonalizationAlgorithm because
    iteration_stop_tolerance is not a very useful stopping criterion.
    NOTE: this stopping criterion may be removed in the future...
* INTERNAL SOFTWARE MODIFICATIONS:
  - Moved LEVEL_SET_FCN_TYPE and SPATIAL_DERIVATIVE_TYPE outside of the 
    LevelSetMethodToolbox class definition to clean up the code.  They
    are now referred to by using the LSMLIB namespace (e.g. LSMLIB::ENO). 
* BUG FIXES:
  - Fixed arguments passed to ReinitializationAlgorithm and 
    OrthogonalizationAlgorithm in the LevelSetFunctionIntegrator class.
  - Removed LDFLAG* variables from Makefiles for examples since they were
    empty and could potentially override LDFLAG* variables in other included
    Makefiles.
  - Removed some C++ style comments and variable declarations in FMM_* code
    that cause errors on some compilers.

-------------------------------------------------------------------------------
v0.7.0 (2006/02/01)
===================
* ADDED FEATURES:
  - Added methods in LevelSetMethodAlgorithm to create FieldExtensionAlgorithm
    objects.
  - Added accessor method for the TVD Runge-Kutta order to the 
    LevelSetMethodAlgorithm and LevelSetFunctionIntegrator classes.
  - Added cfl_number as an input parameter for the FieldExtensionAlgorithm, 
    ReinitializationAlgorithm, and OrthogonalizationAlgorithm.
  - Modified FieldExtensionAlgorithm so that the user can specify
    which component of level set function to use in field extension 
    calculation.
  - Improved design and performance of LevelSetFunctionIntegrator, 
    ReinitializationAlgorithm, and OrthogonalizationAlgorithm  by 
    removing some redundant PatchData and communications.  The extra 
    communication only occurs during vector level set calculations.
  - Added some debugging output to LevelSetFunctionIntegrator class
    which is activated when verbose_mode is set to true.
  - Redesigned LevelSetMethodAlgorithm to limit the requirements on 
    the LevelSetFunctionIntegratorStrategy class.
  - Added BoundaryConditionToolbox class to make it easier to impose
    common boundary conditions on level set functions.  Support is 
    currently provided for periodic/anti-periodic boundary conditions.
  - Added anti-periodic boundary conditions to LevelSetFunctionIntegrator,
    ReinitializationAlgorithm, OrthogonalizationAlgorithm, and
    FieldExtensionAlgorithm so that they work correctly at periodic
    boundaries (i.e. do not introduce artificial zero level set into
    solution).
  - Improved doxygen documentation that groups public methods together
    by functionality
  - Ported to BlueGene/L with XL compilers.
* INTERNAL SOFTWARE MODIFICATIONS
  - Added LSM3D_ prefix to findLineInTetrahedron() function.
    which is activated when verbose_mode is set to true.
  - Added configuration for different versions of HDF5.
  - Added use of explicit qualifiers for namespaces in explicit template
    instantiation to avoid warnings on some compilers.
  - Improved performance of fortran routines by replacing division with
    multiplication by inverse when possible.  This optimization may not be
    useful for all platforms or compiler option settings.
* BUG FIXES:
  - Fixed bugs in LevelSetFunctionIntegrator when using vector level sets.
  - Fixed bug in computation of stable time step when the external velocity
    field depends on the level set functions.
  - Fixed bugs in FieldExtensionAlgorithm, ReinitializationAlgorithm, and
    OrthogonalizationAlgorithm classes when stop tolerance is not
    specified.
  - Fixed bugs in computation of ENO3 and WENO5 spatial derivatives in 3d.
  - Fixed possible floating-point exception in 
    lsm*dComputeStableNormalVelDt() functions caused by division by zero
    when |grad(phi)| is near zero.

-------------------------------------------------------------------------------
