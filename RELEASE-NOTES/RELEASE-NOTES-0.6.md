Release Notes: LSMLIB (0.6)
===========================

-------------------------------------------------------------------------------
v0.6.3 (2006/01/17)
===================
* BUG FIX:  computation of volume and surface integrals now correctly returns
            the sum reduction across processors.

-------------------------------------------------------------------------------
v0.6.2 (2006/01/13)
===================
* Removed unneeded methods from FieldExtensionAlgorithm, 
  ReinitializationAlgorith, and OrthogonalizationAlgorithm.
* Fixed error reporting in classes in the parallel package. 
* Added resetHierarchyConfiguration() method to LevelSetMethodAlgorithm
  class.
* Changed default spatial derivative order to depend on the spatial 
  derivative type for the LevelSetFunctionIntegrator, 
  ReinitializationAlgorithm, OrthogonalizationAlgorithm, and 
  FieldExtensionAlgorithm classes.
* Changed LevelSetFunctionIntegrator::computeStableDt() so that it 
  ignores user-specified dt values when they are zero or negative.
* Fixed several bugs throughout the library.

-------------------------------------------------------------------------------
v0.6.1 (2006/01/09)
===================
* Added 2nd- and 4th-order central approximations for spatial derivatives.
* Added some const qualifiers to the classes in the parallel package.

-------------------------------------------------------------------------------
v0.6.0 (2006/01/02)
===================
* Reorganized C++ interface so that all classes are in the LSMLIB
  namespace.
* Added support for vector level set method computations.
* Added support for extending fields off of the zero level set via 
  the field extension equation (see Osher & Fedkiw).
* Moved the reinitialization and orthogonalization algorithms to their
  own classes (rather than have them embedded within the 
  LevelSetMethodIntegrator class).
* Renamed methods in the LevelSetMethodPatchStrategy to make their 
  functionality more apparent.
* Added several accessor methods to LevelSetMethodAlgorithm, 
  LevelSetMethodIntegratorStrategy, and LevelSetMethodIntegrator:
  - spatial derivative type and order used in a level set method calculation
  - PatchData handle for control volume data
* Added utility functions to compute integrals over the interior and surface
  of the zero level set of the level set function.
* Removed pointer to GridGeometry from the argument list for the
  LevelSetMethodAlgorithm and LevelSetMethodIntegrator classes because
  it is redundant.
* Added component arguments (with default value 0) to several 
  LevelSetMethodToolbox methods so that they can be used on PatchData 
  with multiple components.
* Added support for doxygen source code documentation generation.
* Minor performance enhancements in numerical kernels (e.g. minimizing
  divisions, etc.)
* Minor "safety" improvements to the C++ classes. 

-------------------------------------------------------------------------------
