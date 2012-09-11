from lsmlib import computeDistanceFunction2d_
from lsmlib import computeExtensionFields2d_
from lsmlib import solveEikonalEquation2d_
import numpy as np

def toFloatArray(arr):
    if type(arr) is not np.ndarray: 
        arr = np.array(arr)
    return arr.astype(float)

def getShape(phi, dx):

    phi = toFloatArray(phi)

    shape = phi.shape

    if type(dx) in (int, float, tuple):
        dx = np.array(dx)

    if dx.shape == ():
        dx = np.resize(dx, (len(shape),))

    assert len(dx) == len(shape), '`phi` and `dx` have incompatible shapes'

    dx = dx.astype(float)
    
    if len(shape) == 1:
        nx, ny = shape[0], 1
        dx, dy = dx[0], 1.
    elif len(shape) == 2:
        ny, nx = shape
        dy, dx = dx[0], dx[1]
    else:
        raise Exception, "3D meshes not yet implemented"

    return nx, ny, dx, dy, shape, phi

def computeDistanceFunction(phi, dx=1., order=2):
    r"""

    Return the distance from the zero contour of the array $\phi$ so
    that $\phi$ satisfies,

    .. math::

       \abs{\nabla \phi} = 1

    :Parameters:

      - `phi`: array-like the zero contour of this array is the
        boundary location for the distance calculation. Phi can of
        1,2,3 or higher dimension
      - `dx`: float or an array-like of shape `len(phi)`, optional the
        cell length in each dimension
      - `order`: int, optional order of computational stencil to use
        in updating points during the fast marching method. Must be 1
        or 2, the default is 2
    
    :Returns:

      - `d` : an array the same shape as phi contains the distance
        from the zero contour (zero level set) of phi to each point in
        the array.

    """
    nx, ny, dx, dy, shape, phi = getShape(phi, dx)
    return computeDistanceFunction2d_(phi.flatten(), nx=nx, ny=ny, dx=dx, dy=dy, order=order).reshape(shape)

def computeExtensionFields(phi, extensionFields, mask=None, extensionMask=None, dx=1., order=2):
    r"""

    Calculates multiple extension fields and distance function

    Calculates extension fields $u$ such that,

    .. math::

      \nabla \phi \cdot \nabla u = 0

    :Parameters:

      - `phi`: array-like the zero contour of this array is the
        boundary location for the distance calculation. Phi can of
        1,2,3 or higher dimension
      - `extensionFields`: array-like extension fields. Extension
        fields have shape `(N,) + s` where N is the number of
        extension fields and `s` is the shape of `phi.shape`
      - `mask`:
      - `extensionMask`: array-like mask of type bool with shape of
        `phi.shape`. Only unmasked values are calculated, but masked
        values are used for initialization
      - `dx`: float or an array-like of shape `len(phi)`, optional the
        cell length in each dimension
      - `order`: int, optional order of computational stencil to use
        in updating points during the fast marching method. Must be 1
        or 2, the default is 2
    
    :Returns:

      - `d, ext` : a tuple containing the calculated distance function
        and extension fields

    """

    nx, ny, dx, dy, shape, phi = getShape(phi, dx)
    extensionFields = toFloatArray(extensionFields)
    extshape = extensionFields.shape

    if extensionFields.shape == shape:
        extensionFields = np.reshape(extensionFields, (1,) + extensionFields.shape)

    assert shape == extensionFields[0].shape, '`phi` and `extensionFields` have incompatible shapes'
        
    if mask is None:
        mask = np.empty((0,), dtype=float)
    else:
        mask = toFloatArray(mask)

        assert shape == mask.shape, '`phi` and `mask` have incompatible shapes'

    if extensionMask is None:
        extensionMask = np.empty((0,), dtype=float)
    else:
        extensionMask = toFloatArray(extensionMask)
        extensionMask = -extensionMask * 2. + 1.
        assert shape == extensionMask.shape, '`phi` and `extensionMask` have incompatible shapes'

    numberOfExtensionFields = extensionFields.shape[0]
    phi, extensionFields = computeExtensionFields2d_(phi.flatten(),
                                                     extensionFields.reshape((numberOfExtensionFields, -1)),
                                                     mask=mask.flatten(),
                                                     extension_mask=extensionMask.flatten(),
                                                     nx=nx, ny=ny, dx=dx, dy=dy, order=order)
                                                     
    return phi.reshape(shape), extensionFields.reshape(extshape)

def solveEikonalEquation(phi, speed, dx=1.):
    r"""

    Calculates `phi` based on speed function `F`

    .. math::

      \abs{\nabla \phi} F = 1

    :Parameters:

      - `phi`: array-like the zero contour of this array is the
        boundary location for the eikonal calculation
      - `speed`: array-like speed function with shape `phi.shape`
      - `dx`: float or an array-like of shape `len(phi)`, optional the
        cell length in each dimension
    
    :Returns:

      - `d`: result of solving the eikonal equation based on $F$

    """

    nx, ny, dx, dy, shape, phi = getShape(phi, dx)
    speed = toFloatArray(speed)
    assert shape == speed.shape

    return solveEikonalEquation2d_(phi.flatten(), speed.flatten(), nx=nx, ny=ny, dx=dx, dy=dy).reshape(shape)

def testing():
    r"""
    Here we will define a few test cases. Firstly a 1D test case

    >>> print np.allclose(computeDistanceFunction((-1., -1., -1., -1., 1., 1., 1., 1.), dx=.5),
    ...                   (-1.75, -1.25, -.75, -0.25, 0.25, 0.75, 1.25, 1.75))
    True

    A 1D test case with very small dimensions.

    >>> dx = 1e-10
    >>> print np.allclose(computeDistanceFunction((-1., -1., -1., -1., 1., 1., 1., 1.), dx=dx),
    ...                   np.arange(8) * dx - 3.5 * dx)
    True

    A 2D test case to test trial values for a pathological case.

    >>> dx = 1.
    >>> dy = 2.
    >>> vbl = -dx * dy / np.sqrt(dx**2 + dy**2) / 2.
    >>> vbr = dx / 2
    >>> vml = dy / 2.
    >>> crossProd = dx * dy
    >>> dsq = dx**2 + dy**2
    >>> top = vbr * dx**2 + vml * dy**2
    >>> sqrt = crossProd**2 *(dsq - (vbr - vml)**2)
    >>> sqrt = np.sqrt(max(sqrt, 0))
    >>> vmr = (top + sqrt) / dsq
    >>> print np.allclose(computeDistanceFunction(((-1., 1., -1.), (1., 1., 1.)), dx=(dx, dy), order=1),
    ...                   ((vbl, vml, vbl), (vbr, vmr, vbr)))
    True

    The `extendVariable` method solves the following equation for a given
    extensionVariable.

    .. math::

       \nabla u \cdot \nabla \phi = 0

    using the fast marching method with an initial condition defined at
    the zero level set.

    >>> tmp = 1 / np.sqrt(2)
    >>> phi = np.array([[-1., 1.], [1., 1.]])
    >>> phi, ext =  computeExtensionFields(phi,
    ...                                    [[-1, .5], [2., -1.]],  
    ...                                    extensionMask=phi < 0,
    ...                                    dx=1., order=1)
    >>> print np.allclose(phi, ((-tmp / 2, 0.5), (0.5, 0.5 + tmp)))
    True
    >>> print np.allclose(ext, [[1.25, .5], [2., 1.25]])
    True

    >>> phi = np.array(((-1., 1., 1.), (1., 1., 1.), (1., 1., 1.)))
    >>> phi, ext = computeExtensionFields(phi,
    ...                                   ((-1., 2., -1.),
    ...                                    (.5, -1., -1.),
    ...                                    (-1., -1., -1.)),
    ...                                   extensionMask=phi < 0,
    ...                                   order=1)
    >>> v1 = 0.5 + tmp
    >>> v2 = 1.5
    >>> tmp1 = (v1 + v2) / 2 + np.sqrt(2. - (v1 - v2)**2) / 2
    >>> tmp2 = tmp1 + 1 / np.sqrt(2)
    >>> print np.allclose(phi, ((-tmp / 2, 0.5, 1.5),
    ...                         (0.5, 0.5 + tmp, tmp1),
    ...                         (1.5, tmp1, tmp2)))
    True
    >>> print np.allclose(ext, ((1.25, 2., 2.),
    ...                         (.5, 1.25, 1.5456),
    ...                         (.5, 0.9544, 1.25)),
    ...                   rtol = 1e-4)
    True

    Test case for a bug that occurs when initializing the distance
    variable at the interface. Currently it is assumed that adjacent
    cells that are opposite sign neighbors have perpendicular normal
    vectors. In fact the two closest cells could have opposite
    normals.

    >>> print np.allclose(computeDistanceFunction((-1., 1., -1.)), (-0.5, 0.5, -0.5))
    True

    Testing second order. This example failed with Scikit-fmm.

    >>> phi = ((-1., -1., 1., 1.),
    ...        (-1., -1., 1., 1.),
    ...        (1., 1., 1., 1.),
    ...        (1., 1., 1., 1.))
    >>> answer = ((-1.30473785, -0.5, 0.5, 1.49923009),
    ...           (-0.5, -0.35355339, 0.5, 1.45118446),
    ...           (0.5, 0.5, 0.97140452, 1.76215286),
    ...           (1.49923009, 1.45118446, 1.76215286, 2.33721352))
    >>> print np.allclose(computeDistanceFunction(phi),
    ...                   answer,
    ...                   rtol=1e-9)
    True

    Solve the level set equation in two dimensions for a circle. 

    The 2D level set equation can be written,

    .. math::

        \abs{\nabla \phi} = 1

    and the boundary condition for a circle is given by, :math:`\phi = 0` at 
    :math:`(x - L / 2)^2 + (y - L / 2)^2 = (L / 4)^2`.

    The solution to this problem will be demonstrated in the following
    script. Firstly, setup the parameters.

    >>> def mesh(nx=1, ny=1, dx=1., dy=1.):
    ...     y, x = np.mgrid[0:nx,0:ny]
    ...     x = x * dx + dx / 2
    ...     y = y * dy + dy / 2
    ...     return x, y

    >>> dx = 1.
    >>> N = 11
    >>> L = N * dx
    >>> x, y = mesh(nx=N, ny=N, dx=dx, dy=dx) 
    >>> phi = -np.ones(N * N, 'd')
    >>> phi[(x.flatten() - L / 2.)**2 + (y.flatten() - L / 2.)**2 < (L / 4.)**2] = 1.
    >>> phi = np.reshape(phi, (N, N))
    >>> phi = computeDistanceFunction(phi, dx=dx, order=1).flatten()

    >>> dX = dx / 2.
    >>> m1 = dX * dX / np.sqrt(dX**2 + dX**2)
    >>> def evalCell(phix, phiy, dx):
    ...     aa = dx**2 + dx**2
    ...     bb = -2 * ( phix * dx**2 + phiy * dx**2)
    ...     cc = dx**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dx**2
    ...     sqr = np.sqrt(bb**2 - 4. * aa * cc)
    ...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
    >>> v1 = evalCell(-dX, -m1, dx)[0] 
    >>> v2 = evalCell(-m1, -dX, dx)[0]
    >>> v3 = evalCell(m1,  m1,  dx)[1]
    >>> v4 = evalCell(v3, dX, dx)[1]
    >>> v5 = evalCell(dX, v3, dx)[1]
    >>> MASK = -1000.
    >>> trialValues = np.array((   
    ...     MASK,  MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK,
    ...     MASK,  MASK, MASK, MASK,-3*dX,-3*dX,-3*dX, MASK, MASK, MASK, MASK,
    ...     MASK,  MASK, MASK,   v1,  -dX,  -dX,  -dX,   v1, MASK, MASK, MASK,
    ...     MASK,  MASK,   v2,  -m1,   m1,   dX,   m1,  -m1,   v2, MASK, MASK,
    ...     MASK, -dX*3,  -dX,   m1,   v3,   v4,   v3,   m1,  -dX,-dX*3, MASK,
    ...     MASK, -dX*3,  -dX,   dX,   v5, MASK,   v5,   dX,  -dX,-dX*3, MASK,
    ...     MASK, -dX*3,  -dX,   m1,   v3,   v4,   v3,   m1,  -dX,-dX*3, MASK,
    ...     MASK,  MASK,   v2,  -m1,   m1,   dX,   m1,  -m1,   v2, MASK, MASK,
    ...     MASK,  MASK, MASK,   v1,  -dX,  -dX,  -dX,   v1, MASK, MASK, MASK,
    ...     MASK,  MASK, MASK, MASK,-3*dX,-3*dX,-3*dX, MASK, MASK, MASK, MASK,
    ...     MASK,  MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK), 'd')

    >>> phi[trialValues == MASK] = MASK
    >>> print np.allclose(phi, trialValues)
    True

    Here we solve the level set equation in two dimensions for a square. The equation is
    given by:

    .. math::

       \abs{\nabla \phi} &= 1 \\
       \phi &= 0 \qquad \text{at} \qquad \begin{cases}
           x = \left( L / 3, 2 L / 3 \right) 
           & \text{for $L / 3 \le y \le 2 L / 3$} \\
           y = \left( L / 3, 2 L / 3 \right) 
           & \text{for $L / 3 \le x \le 2 L / 3$}
       \end{cases}

    Do the tests:

    >>> dx = 0.5
    >>> dy = 2.
    >>> nx = 5
    >>> ny = 5
    >>> Lx = nx * dx
    >>> Ly = ny * dy

    >>> x, y = mesh(nx=nx, ny=ny, dx=dx, dy=dy)
    >>> x = x.flatten()
    >>> y = y.flatten()
    >>> phi = -np.ones(nx * ny, 'd')
    >>> phi[((Lx / 3. < x) & (x < 2. * Lx / 3.)) & ((Ly / 3. < y) & (y < 2. * Ly / 3))] = 1.
    >>> phi = np.reshape(phi, (nx, ny))
    >>> phi = computeDistanceFunction(phi, dx=(dy, dx), order=1).flatten()

    >>> def evalCell(phix, phiy, dx, dy):
    ...     aa = dy**2 + dx**2
    ...     bb = -2 * ( phix * dy**2 + phiy * dx**2)
    ...     cc = dy**2 * phix**2 + dx**2 * phiy**2 - dx**2 * dy**2
    ...     sqr = np.sqrt(bb**2 - 4. * aa * cc)
    ...     return ((-bb - sqr) / 2. / aa,  (-bb + sqr) / 2. / aa)
    >>> val = evalCell(-dy / 2., -dx / 2., dx, dy)[0]
    >>> v1 = evalCell(val, -3. * dx / 2., dx, dy)[0]
    >>> v2 = evalCell(-3. * dy / 2., val, dx, dy)[0]
    >>> v3 = evalCell(v2, v1, dx, dy)[0]
    >>> v4 = dx * dy / np.sqrt(dx**2 + dy**2) / 2
    >>> arr = np.array((
    ...     v3           , v2      , -3. * dy / 2.   , v2      , v3,
    ...     v1           , val     , -dy / 2.        , val     , v1           ,
    ...     -3. * dx / 2., -dx / 2., v4              , -dx / 2., -3. * dx / 2.,
    ...     v1           , val     , -dy / 2.        , val     , v1           ,
    ...     v3           , v2      , -3. * dy / 2.   , v2      , v3           ))
    >>> print np.allclose(arr, phi)
    True

    Check some assertion errors

    >>> computeDistanceFunction([[-1, 1],[1, 1]], dx=(1, 2, 3))
    Traceback (most recent call last):
      ...
    AssertionError: `phi` and `dx` have incompatible shapes
    >>> computeExtensionFields([[-1, 1],[1, 1]], extensionFields=[1, 1])
    Traceback (most recent call last):
      ...
    AssertionError: `phi` and `extensionFields` have incompatible shapes
    """

    pass

distance = computeDistanceFunction

def test():
    r"""
    Run all the doctests available.
    """
    import doctest
    import pylsmlib
    doctest.testmod(pylsmlib)
