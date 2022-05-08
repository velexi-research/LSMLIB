Release Notes: LSMLIB v0.5
==========================

-------------------------------------------------------------------------------
v0.5.3 (2005/10/25)
===================
* Redesigned the LevelSetMethod classes to make level set function data 
  available to the LevelSetMethodVelocityFieldStrategy class when 
  computing and initializing the velocity field.
* Added 2d example for LevelSetMethod classes when simulating motion under 
  a normal velocity.
* Minor bug fixes and corrections in comments/documentation.
* Split the samrai examples into two separate types: advection and normal 
  velocity motion.
* Added mex files to matlab package
  - 3D upwind-HJ spatial derivative calculations
  - 3D fast marching method calculations

-------------------------------------------------------------------------------
v0.5.2 (2005/10/20)
===================
* Added support for motion resulting from specification of the normal 
  velocity.
  - This changed the LevelSetMethodIntegrator, 
    LevelSetMethodVelocityFieldStrategy and LevelSetMethodToolBox classes.
* Added 3d example for LevelSetMethod classes.

-------------------------------------------------------------------------------
v0.5.1 (2005/10/12)
===================
* Added restart capabilities for the LevelSetMethod classes.
* Added getIntegratorStep() method for LevelSetMethodAlgorithm, 
  LevelSetMethodIntegratorStrategy, and LevelSetMethodIntegrator
  classes to keep track of how many integration time steps have been 
  taken.
* Added setCurrentTime() method to the LevelSetMethodVelocityFieldStrategy
  to make it possible to synchronize the simulation time of the velocity 
  field strategy with the simulation time of the LevelSetMethodIntegrator.
* Fixed various bugs in the code (primarily parallel bugs).

-------------------------------------------------------------------------------
v0.5.0 (2005/10/11)
===================
* Initial full featured version of parallel, level set method library 
  for codimension-one problems.

-------------------------------------------------------------------------------
