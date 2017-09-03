Gas Bubble Simulation
=====================

This example uses the parallel LSMLIB package to simulate gas bubble evolution.
The normal velocity at the interface is given by

  ```
  (K * pi / V - P_ext) - sigma * kappa
  ```

where kappa is the mean curvature of the interface, V is the volume of the
bubble, and K, P_ext and sigma are constants respectively representing the
relationship between the pressure of the gas contained within the bubble
(modeled as an ideal gas), P_ext is the external ambient pressure, and sigma
is the surface tension.

The level see simulation for this example utilizes the following LSMLIB
capabilities:

- numerical computation of mean curvature term
- level set evolution using a normal velocity term
- extension of the normal velocity off of the interface
- periodic reinitialization of level set function
- computation of volume defined by level set function 

Accuracy of the simulation is controlled by setting the method and order
of the spatial derivatives and time integration.

Model Notes
-----------
This physical model used by this model of interface motion driven by
a combination pressure and surface tension forces.

Building Example
----------------
1. Build LSMLIB library (see instructions in INSTALL.markdown).

2. Build example.

   ```shell
   $ cd ${BUILD_DIR}
   $ make gas-bubble-parallel
   ```
Running Example
---------------
1. Run executable.

   * Serial

   ```shell
   $ cd ${EXAMPLE_DIR}
   $ ./gas-bubble-parallel <INPUT_FILE>
   ```

   * Parallel

   ```shell
   $ cd ${EXAMPLE_DIR}
   $ mpirun -np <NUM_PROCS> gas-bubble-parallel <INPUT_FILE>
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
