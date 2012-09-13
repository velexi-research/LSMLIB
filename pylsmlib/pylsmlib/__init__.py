from lsmlib import computeDistanceFunction2d_
from lsmlib import computeExtensionFields2d_
from lsmlib import solveEikonalEquation2d_
import numpy as np

__docformat__ = 'restructuredtext'

def _getVersion():
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        version = get_distribution(__name__).version
    except DistributionNotFound:
        version = "unknown, try running `python setup.py egg_info`"
        
    return version

__version__ = _getVersion()

def toFloatArray(arr):
    if type(arr) is not np.ndarray: 
        arr = np.array(arr)
    return arr.astype(float)

def getShape(phi, dx, order):

    phi = toFloatArray(phi)

    shape = phi.shape

    if type(dx) in (int, float, tuple):
        dx = np.array(dx)

    if dx.shape == ():
        dx = np.resize(dx, (len(shape),))

    assert len(dx) == len(shape), '`phi` and `dx` have incompatible shapes'

    assert order in (1, 2), '`order` is 1 or 2'

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

def computeDistanceFunction(phi0, dx=1., order=2):
    r"""

    Solves

    .. math::

       |\nabla \phi| = 1 \;\; \text{with} \;\; \phi=0 \;\; \text{where} \;\; \phi_0=0

    for :math:`\phi` using the zero level set of :math:`\phi_0`.
    
    :Parameters:

      - `phi0`: array of positive and negative values locating the
        zero level set of :math:`\phi`, shape determines the
        dimension.
      - `dx`: the cell dimensions
      - `order`: order of the computational stencil, either 1 or 2
    
    :Returns:

      - the calculated distance function, :math:`\phi`

    """
    nx, ny, dx, dy, shape, phi0 = getShape(phi0, dx, order)
    return computeDistanceFunction2d_(phi0.flatten(), nx=nx, ny=ny, dx=dx, dy=dy, order=order).reshape(shape)

def computeExtensionFields(phi0, u0, mask=None, u0mask=None, dx=1., order=2):
    r"""

    Solves

    .. math::

      \nabla \phi \cdot \nabla u = 0 \;\; \text{with} \;\; u=u_0 \;\; \text{where} \;\; \phi_0=0 

    for multiple extension fields :math:`u` and
  
    .. math::

      |\nabla \phi| = 1 \;\; \text{with} \;\; \phi=0 \;\; \text{where} \;\; \phi_0=0

    for :math:`\phi` using the zero level set of :math:`\phi_0`.

    :Parameters:


      - `phi0`: array of positive and negative values locating the
        zero level set of :math:`\phi`, shape determines the
        dimension.
      - `u0`: array of values used as the initial condition to
        calculate the extension fields, :math:`u`; shape is either
        `phi0.shape` or `(N,) + phi0.shape`, where `N` is the number
        of fields to extend
      - `mask`:
      - `u0mask`: array of Boolean values determining which values of
        :math:`u` are calculated; `True` values set :math:`u=u_0`,
        `False` values set :math:`u=u` determined from the eikonal
        equation; shape is `phi0.shape`
      - `dx`: the cell dimensions
      - `order`: order of the computational stencil, either 1 or 2
    
    :Returns:

      - a tuple containing (:math:`\phi`, :math:`u`)

    """

    nx, ny, dx, dy, shape, phi0 = getShape(phi0, dx, order)
    u0 = toFloatArray(u0)
    u0shape = u0.shape

    if u0shape == shape:
        u0 = np.reshape(u0, (1,) + u0shape)

    assert shape == u0[0].shape, '`phi` and `extensionFields` have incompatible shapes'
        
    if mask is None:
        mask = np.empty((0,), dtype=float)
    else:
        mask = toFloatArray(mask)

        assert shape == mask.shape, '`phi` and `mask` have incompatible shapes'

    if u0mask is None:
        u0mask = np.empty((0,), dtype=float)
    else:
        u0mask = toFloatArray(u0mask)
        u0mask = -u0mask * 2. + 1.
        assert shape == u0mask.shape, '`phi` and `extensionMask` have incompatible shapes'

    N = u0.shape[0]
    phi, u = computeExtensionFields2d_(phi0.flatten(),
                                       u0.reshape((N, -1)),
                                       mask=mask.flatten(),
                                       extension_mask=u0mask.flatten(),
                                       nx=nx, ny=ny, dx=dx, dy=dy, order=order)
                                                     
    return phi.reshape(shape), u.reshape(u0shape)

def solveEikonalEquation(phi0, speed, dx=1.):
    r"""

    Calculates :math:`\phi` based on the speed function, :math:`F`,
    using

    .. math::

      |\nabla \phi| F = 1 \;\; \text{with} \;\; \phi=0 \;\; \text{where} \;\; \phi_0=0

    :Parameters:

      - `phi0`: array of positive and negative values locating the
        zero level set of :math:`\phi`, shape determines the
        dimension.
      - `speed`: speed function, :math:`F`, shape is `phi0.shape`
      - `dx`: the cell dimensions
    
    :Returns:

      - the field calculated by solving the eikonal equation,
        :math:`\phi`

    """

    nx, ny, dx, dy, shape, phi0 = getShape(phi0, dx, 1)
    speed = toFloatArray(speed)
    assert shape == speed.shape

    return solveEikonalEquation2d_(phi0.flatten(), speed.flatten(), nx=nx, ny=ny, dx=dx, dy=dy).reshape(shape)

def testing():
    r"""

    **1D Test**

    >>> print np.allclose(computeDistanceFunction((-1., -1., -1., -1., 1., 1., 1., 1.), dx=.5),
    ...                   (-1.75, -1.25, -.75, -0.25, 0.25, 0.75, 1.25, 1.75))
    True

    Small dimensions.

    >>> dx = 1e-10
    >>> print np.allclose(computeDistanceFunction((-1., -1., -1., -1., 1., 1., 1., 1.), dx=dx),
    ...                   np.arange(8) * dx - 3.5 * dx)
    True

    **Bug Fix**
  
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

    **Test Extension Field Calculation**

    >>> tmp = 1 / np.sqrt(2)
    >>> phi = np.array([[-1., 1.], [1., 1.]])
    >>> phi, ext =  computeExtensionFields(phi,
    ...                                    [[-1, .5], [2., -1.]],  
    ...                                    u0mask=phi < 0,
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
    ...                                   u0mask=phi < 0,
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

    **Bug Fix**
 
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

    **Circle Example** 

    Solve the level set equation in two dimensions for a circle. 

    The 2D level set equation can be written,

    .. math::

        |\nabla \phi| = 1

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

    **Square Example**

    Here we solve the level set equation in two dimensions for a square. The equation is
    given by:

    .. math::

       |\nabla \phi| &= 1 \\
       \phi &= 0 \qquad \text{at} \qquad \begin{cases}
           x = \left( L / 3, 2 L / 3 \right) 
           & \text{for $L / 3 \le y \le 2 L / 3$} \\
           y = \left( L / 3, 2 L / 3 \right) 
           & \text{for $L / 3 \le x \le 2 L / 3$}
       \end{cases}

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

    **Assertion Errors**

    >>> computeDistanceFunction([[-1, 1],[1, 1]], dx=(1, 2, 3))
    Traceback (most recent call last):
      ...
    AssertionError: `phi` and `dx` have incompatible shapes
    >>> computeExtensionFields([[-1, 1],[1, 1]], u0=[1, 1])
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
