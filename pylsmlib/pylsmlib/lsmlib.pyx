import cython
import numpy as np
from numpy.core import intc
from numpy.compat import asbytes
from numpy import int, double, int32
from lsmlib cimport computeDistanceFunction3d
from lsmlib cimport computeExtensionFields3d
from lsmlib cimport solveEikonalEquation3d
from libc.stdlib cimport malloc, free

def computeDistanceFunction3d_(np.ndarray[double, ndim=1] phi, nx=1, ny=1, nz=1, dx=1., dy=1., dz=1., order=2):
    cdef np.ndarray[double, ndim=1] distance_function = np.zeros((nx * ny * nz,))
    cdef int spatial_derivative_order = order
    cdef np.ndarray[int, ndim=1] grid_dims = np.array((nx, ny, nz), dtype=int32)
    cdef np.ndarray[double, ndim=1] _dx = np.array((dx, dy, dz))

    cdef double *maskdata=NULL

    error = computeDistanceFunction3d(
		<double *> distance_function.data,
    	<double *> phi.data,
    	<double *> maskdata,
    	spatial_derivative_order,
    	<int *> grid_dims.data,
    	<double *> _dx.data)

    return distance_function

def computeExtensionFields3d_(np.ndarray[double, ndim=1] phi, np.ndarray[double, ndim=2] extensionFields, np.ndarray[double, ndim=1] mask=np.empty((0,), dtype=double), np.ndarray[double, ndim=1] extension_mask=np.empty((0,), dtype=double), nx=1,  ny=1, nz=1, dx=1., dy=1., dz=1., order=2):

    cdef int num_ext_fields = extensionFields.shape[0]

    cdef np.ndarray[double, ndim=1] distance_function = np.zeros((nx * ny * nz,))
    cdef int spatial_derivative_order = order

    cdef np.ndarray[int, ndim=1] grid_dims = np.array((nx, ny, nz), dtype=int32)
    cdef np.ndarray[double, ndim=1] _dx = np.array((dx, dy, dz))

    cdef double **ext_fields = <double **> malloc(num_ext_fields * sizeof(double*))
    cdef double **source_fields = <double **> malloc(num_ext_fields * sizeof(double*))
    cdef np.ndarray[double, ndim=2] extReturnFields = np.zeros((num_ext_fields, nx * ny * nz))

    cdef double *maskdata=NULL
    cdef double *extension_maskdata=NULL

    if len(mask) > 0:
        maskdata = <double *> mask.data

    if len(extension_mask) > 0:
        extension_maskdata = <double *> extension_mask.data


    for i in range(num_ext_fields):
        ext_fields[i] = &extReturnFields[i,0]
        source_fields[i] = &extensionFields[i,0]

    error = computeExtensionFields3d(
		<double *> distance_function.data,
		<double **> ext_fields,
    	<double *> phi.data,
    	<double *> maskdata,
		<double **> source_fields,
		<double *> extension_maskdata,
		num_ext_fields,
    	spatial_derivative_order,
    	<int *> grid_dims.data,
    	<double *> _dx.data)

    free(ext_fields)
    free(source_fields)

    return distance_function, extReturnFields

def solveEikonalEquation3d_(np.ndarray[double, ndim=1] phi,
                            np.ndarray[double, ndim=1] speed, nx=1, ny=1, nz=1, dx=1., dy=1., dz=1., order=2):

    cdef np.ndarray[double, ndim=1] mask = np.zeros((nx * ny * nz,))
    cdef int spatial_derivative_order = order
    cdef np.ndarray[int, ndim=1] grid_dims = np.array((nx, ny, nz), dtype=int32)
    cdef np.ndarray[double, ndim=1] _dx = np.array((dx, dy, dz))

    cdef double *maskdata=NULL

    minphi = min(phi)
    phi[phi > 0] = minphi - 1.
    phi[:] = phi - minphi

    error = solveEikonalEquation3d(
		<double *> phi.data,
    	<double *> speed.data,
    	<double *> maskdata,
    	spatial_derivative_order,
    	<int *> grid_dims.data,
    	<double *> _dx.data)

    phi[:] = phi + minphi

    return phi
