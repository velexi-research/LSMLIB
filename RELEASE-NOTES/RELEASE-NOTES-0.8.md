Release Notes: LSMLIB (0.8)
===========================

-------------------------------------------------------------------------------
v0.8.0 (2006/06/02)
===================
* ENHANCEMENTS:
  - New support for serial calculations has been added thanks to generous
    contributions by Dr. Masa Prodanovic. 
  - Modified fast marching method algorithm to set ghostcell values by
    taking the value of the nearest interior grid point rather than
    using a default value of zero.  This change was made to avoid the
    introduction of spurious zero level sets.
  - Changed logic for computation of a stable time step to improve performance.
  - MATLAB SUPPORT:  
    * added support for computation of plus/minus HJ ENO/WENO spatial 
      derivatives, TVD Runge-Kutta time integration, advancing level set 
      functions in time using the level set evolution equation.
    * added MEX-file support for Linux, Windows, and Solaris platforms.
  - SAMRAI is no longer a requirement for LSMLIB.  The serial and MATLAB 
    capabilities of LSMLIB are still available when LSMLIB is not configured
    with SAMRAI.
  - The FieldExtensionAlgorithm now supports anti-periodic boundary 
    conditions for the field variable that is being extended off of the
    zero level set.
  - The user's guide has been expanded (but will continue to grow in 
    the future releases).
  - The configure script now builds object code that is position independent
    when the architecture supports it.  This makes it possible to link
    LSMLIB into shared objects.
* API CHANGES:
  - The reinitializeLevelSetFunctions() method for the LevelSetMethodAlgorithm,
    LevelSetFunctionIntegratorStrategy, and LevelSetFunctionIntegrator classes
    has been changed to take an optional argument that specifies which 
    level set function (i.e. PHI or PSI) to reinitialize.  By default, 
    PHI is reinitialized, so codimension-one problems will not need to 
    pass an argument.
  - The FieldExtensionAlgorithm constructors now have an optional argument
    that may be used to specify that the extension field should use 
    anti-periodic boundary conditions.
  - The LevelSetMethodAlgorithm, ReinitializationAlgorithm,
    OrthogonalizationAlgorithm, FieldExtensionAlgorithm, 
    LevelSetFunctionIntegratorStrategy, and LevelSetFunctionIntegrator
    have all been modified so that reinitialization, orthogonalization,
    and computation of the extension field all take an optional max_iterations
    argument that can be used to override the stopping criterion specified
    in the input file.
  - The LevelSetMethodAlgorithm, LevelSetFunctionIntegrator, and 
    LevelSetFunctionIntegratorStrategy classes have been expanded to 
    include accessor methods for the reinitialization and orthogonalization
    interval.  These methods are provided to allow software built using
    LSMLIB to override values specified via the input file.
  - The fast marching method toolbox functions have been modified so that
    it is NOT necessary to pass ghostbox/ghostcell width information.  At
    grid points that border the edge of the computational domain, the values 
    of the distance function and extension fields are computed using only 
    the values at interior grid points.  This choice could potentially lead
    to some errors near the boundaries of the computational domain but should 
    not be a big issue for problems where the zero level set is sufficiently
    far from the boundaries. 
* MINOR IMPROVEMENTS:
  - Improved AntiPeriodicBoundaryConditionModule to be more robust about 
    detecting anti-periodic boundary conditions.
  - Improved the autoconf set_build_mode.m4 macro.
  - Added extern "C" guards around C header files.
  - Improved pattern rules in Makefile.config.in to avoid conflicts with
    compiler flags set by other software libraries.
* BUG FIXES/INTERNAL SOFTWARE MODIFICATIONS:
  - Fixed bug in OrthogonalizationAlgorithm where the wrong data is used
    to compute the gradient of fixed field during the orthogonalization
    procedure.
  - Fixed bugs in numerical kernels for field extension so that they no 
    longer allow grid cells immediately adjacent (along the coordinate
    axes) to the zero level set to be modified during the field extension 
    algorithm.
  - Redesigned OrthogonalizationAlgorithm to internally use 
    FieldExtensionAlgorithm instead of having redundant code.
  - Fixed bug in fast marching method code that caused the front to be 
    incorrectly initialized when grid points lie exactly on the zero level
    set.
  - Reorganized the logic for orthogonalization and reinitialization of
    level set functions.
  - Updated MATLAB MEX files to reflect changes in fast marching method
    and spatial derivative function signatures.
  - Fixed a minor pointer to integer conversion bug in FMM_Heap that
    arises on 64-bit architectures.
  - Fixed some bugs in the configure script.

-------------------------------------------------------------------------------
