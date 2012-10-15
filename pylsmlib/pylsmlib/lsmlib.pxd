import numpy as np
cimport numpy as np

cdef extern int computeDistanceFunction3d(
	double *distance_function,
	double *phi,
	double *mask,
	int spatial_derivative_order,
	int *grid_dims,
	double *dx)

cdef extern int computeExtensionFields3d(
     double *distance_function,
     double **extension_fields,
     double *phi,
     double *mask,
     double **source_fields,
     double *extension_mask,
     int num_ext_fields,
     int spatial_derivative_order,
     int *grid_dims,
     double *dx)

cdef extern int solveEikonalEquation3d(
	double *phi,
	double *speed,
	double *mask,
	int spatial_derivative_order,
	int *grid_dims,
	double *dx)


