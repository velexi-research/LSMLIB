import cython
import numpy as np
from numpy.core import intc	
cimport numpy as np
from numpy.compat import asbytes
from numpy import int, double
from pylsmlib cimport computeDistanceFunction2d

def computeDistanceFunction2d_(phi, nx=1, ny=1, dx=1., dy=1.):
    cdef np.ndarray[double, ndim=1] _phi = phi.copy()
    cdef np.ndarray[double, ndim=1] distance_function = np.zeros((nx * ny,))
    cdef np.ndarray[double, ndim=0] mask = np.array(0.0)
    cdef int spatial_derivative_order = 2							
    cdef np.ndarray[long int, ndim=1] grid_dims = np.array((nx, ny), dtype=int)
    cdef np.ndarray[double, ndim=1] dX = np.array((dx, dy))
    
    error = computeDistanceFunction2d(
	<double *> distance_function.data,
    	<double *> _phi.data,	
    	<double *> mask.data,
    	spatial_derivative_order,
    	<int *> grid_dims.data,
    	<double *> dX.data)

    return distance_function