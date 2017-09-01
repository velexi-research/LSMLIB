Curvature-Driven Interface Motion
=================================

This example uses the parallel LSMLIB package to simulate curvature-driven
interface motion. The normal velocity at the interface is given by

  a - b * kappa

where a and b are constants and kappa is the mean curvature of the interface.
The level set simulation for this example utilizes the following LSMLIB
capabilities:

- numerical computation of mean curvature term
- level set evolution using a normal velocity term
- extension of the normal velocity off of the interface
- periodic reinitialization of level set function

Accuracy of the simulation is controlled by setting the method and order
of the spatial derivatives and time integration.

Model Notes
-----------
This mathematical model is an idealization of interface motion driven by
a combination pressure and surface tension forces. It should be noted, however,
that it _neglects_ critical physical effects that are important for simulating
real physical systems. For example, it _cannot_ be used to study equilibrium
interfaces.

Building Example
----------------
1. Build LSMLIB library (see instructions in INSTALL.markdown).

2. Build example.

   ```shell
   $ cd ${BUILD_DIR}
   $ make curvature
   ```
Running Example
---------------
1. Run executable.

   * Serial

   ```shell
   $ cd ${EXAMPLE_DIR}
   $ ./curvature <INPUT_FILE>
   ```

   * Parallel

   ```shell
   $ cd ${EXAMPLE_DIR}
   $ mpirun -np <NUM_PROCS> curvature <INPUT_FILE>
   ```

Output
------
- CURVATURE_*.restart: restart files

- CURVATURE_*.visit: VisIt visualization files

- CURVATURE_*.log.MPI_NODE_ID: log files

Visualization
-------------
- Use the VisIt visualization tool developed by Lawrence Livermore National
  Laboratory (LLNL) to open the dumps.visit file.
