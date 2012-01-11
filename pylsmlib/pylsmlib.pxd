import numpy as np
cimport numpy as np

cdef extern int computeDistanceFunction2d(
     double *distance_function,
     double *_phi,
     double *mask,
     int spatial_derivative_order,
     int *grid_dims,
     double *dX)