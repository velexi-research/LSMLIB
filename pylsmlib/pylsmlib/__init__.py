from lsmlib import computeDistanceFunction2d_
from lsmlib import computeExtensionFields2d_
from lsmlib import solveEikonalEquation2d_
import numpy as np

def getShape(shape, dx):

    if type(dx) in (int, float, tuple):
        dx = np.array(dx)

    if dx.shape == ():
        dx = np.resize(dx, shape)

    if not len(dx) == len(shape):
        raise Exception, "Arrays are not the correct shape"

    dx = dx.astype(float)
    
    if len(shape) == 1:
        nx, ny = shape[0], 1
        dx, dy = dx[0], 1.
    elif len(shape) == 2:
        ny, nx = shape
        dy, dx = dx[0], dx[1]
    else:
        raise Exception, "3D meshes not yet implemented"

    return nx, ny, dx, dy

def computeDistanceFunction(phi, dx=1., order=2):
    phi = phi.astype(float)
    shape = phi.shape
    nx, ny, dx, dy = getShape(shape, dx)
    return computeDistanceFunction2d_(phi.flatten(), nx=nx, ny=ny, dx=dx, dy=dy, order=order).reshape(shape)

def computeExtensionFields(phi, extensionFields, mask=None, extensionMask=None, dx=1., order=2):
    
    phi = phi.astype(float)
    extensionFields = extensionFields.astype(float)
    mask = mask.astype(float)
    
    if not (phi.shape == extensionFields[0].shape and phi.shape == mask.shape):
        raise Exception, "Arrays are not the correct shape"

    if mask is None:
        mask = np.empty((0,), dtype=float)
    else:
        if not (phi.shape == mask.shape):
            raise Exception, "Arrays are not the correct shape"

    if extensionMask is None:
        extensionMask = np.empty((0,), dtype=float)
    else:
        if not (phi.shape == extensionMask.shape):
            raise Exception, "Arrays are not the correct shape"

    shape = phi.shape
    nx, ny, dx, dy = getShape(shape, dx)
    numberOfExtensionFields = extensionFields.shape[0]
    
    phi, extensionFields = computeExtensionFields2d_(phi.flatten(),
                                                     extensionFields.reshape((numberOfExtensionFields, -1)),
                                                     mask=mask,
                                                     extension_mask=extensionMask,
                                                     nx=nx, ny=ny, dx=dx, dy=dy, order=order)
                                                     
    return phi.reshape(shape), extensionFields.reshape((numberOfExtensionFields,) + shape)

def solveEikonalEquation(phi, speed, dx=1.):
    
    phi = phi.astype(float)
    shape = phi.shape
    nx, ny, dx, dy = getShape(shape, dx)
    
    if not (phi.shape == speed.shape):
        raise Exception, "Arrays are not the correct shape"

    return solveEikonalEquation2d_(phi.flatten(), speed.flatten(), nx=nx, ny=ny, dx=dx, dy=dy).reshape(shape)

distance = computeDistanceFunction
