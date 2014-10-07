from lsmlib import computeDistanceFunction3d_
from lsmlib import computeExtensionFields3d_
from lsmlib import solveEikonalEquation3d_
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

    if type(dx) in (int, float, tuple, list):
        dx = np.array(dx)

    if dx.shape == ():
        dx = np.resize(dx, (len(shape),))

    if len(dx) != len(shape):
        raise ValueError, "dx must be of length len(phi.shape)"

    if order not in (1, 2):
        raise ValueError, 'order must be 1 or 2'

    dx = dx.astype(float)

    if np.prod(dx) == 0:
        raise ValueError, 'dx must be greater than zero'
    
    if len(shape) == 1:
        nx, ny = shape[0], 1
        dx, dy = dx[0], 1.
    elif len(shape) == 2:
        ny, nx = shape
        dy, dx = dx[0], dx[1]
    elif len(shape) == 3:
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
    return computeDistanceFunction3d_(phi0.flatten(), nx=nx, ny=ny, dx=dx, dy=dy, order=order).reshape(shape)

def computeExtensionFields(phi0, speed, dx=1., order=2, mask=None, ext_mask=None, ):
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
      - `speed`: array of values used as the initial condition to
        calculate the extension fields, :math:`u0`; shape is either
        `phi0.shape` or `(N,) + phi0.shape`, where `N` is the number
        of fields to extend
      - `mask`:
      - `ext_mask`: array of Boolean values determining which values of
        :math:`u` are calculated; `True` values set :math:`u=u_0`,
        `False` values set :math:`u=u` determined from the eikonal
        equation; shape is `phi0.shape`
      - `dx`: the cell dimensions
      - `order`: order of the computational stencil, either 1 or 2

    :Returns:

      - a tuple containing (:math:`\phi`, :math:`u`)

    """

    nx, ny, dx, dy, shape, phi0 = getShape(phi0, dx, order)
    u0 = toFloatArray(speed)
    u0shape = u0.shape

    if u0shape == shape:
        u0 = np.reshape(u0, (1,) + u0shape)

    if shape != u0[0].shape:
        raise ValueError, "phi and speed must have the same shape"

    if mask is None:
        mask = np.empty((0,), dtype=float)
    else:
        mask = toFloatArray(mask)

        assert shape == mask.shape, '`phi` and `mask` have incompatible shapes'

    if ext_mask is None:
        ext_mask = np.empty((0,), dtype=float)
    else:
        ext_mask = toFloatArray(ext_mask)
        ext_mask = -ext_mask * 2. + 1.
        assert shape == ext_mask.shape, '`phi` and `extensionMask` have incompatible shapes'

    N = u0.shape[0]
    phi, u = computeExtensionFields3d_(phi0.flatten(),
                                       u0.reshape((N, -1)),
                                       mask=mask.flatten(),
                                       extension_mask=ext_mask.flatten(),
                                       nx=nx, ny=ny, dx=dx, dy=dy, order=order)

    return phi.reshape(shape), u.reshape(u0shape)

def solveEikonalEquation(phi0, speed, dx=1., order=2):
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
      - `order`: order of the computational stencil, either 1 or 2
      
    :Returns:

      - the field calculated by solving the eikonal equation,
        :math:`\phi`

    """

    nx, ny, dx, dy, shape, phi0 = getShape(phi0, dx, 1)
    speed = toFloatArray(speed)

    if shape != speed.shape:
        raise ValueError, "phi and speed must have the same shape"

    return solveEikonalEquation3d_(phi0.flatten(), speed.flatten(), nx=nx, ny=ny, dx=dx, dy=dy, order=order).reshape(shape)

def travel_time(phi0, speed, dx=1., order=2):
    return abs(solveEikonalEquation(phi0, speed, dx, order))

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
    ...                                    ext_mask=phi < 0,
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
    ...                                   ext_mask=phi < 0,
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

    **A test for a bug in both LSMLIB and Scikit-fmm**

    The following test gives different results depending on whether
    LSMLIB_ or Scikit-fmm_ is used. This issue occurs when calculating
    second order accurate distance functions. When a value becomes
    "known" after previously being a "trial" value it updates its
    neighbors' values. In a second order scheme the neighbors one step
    away also need to be updated (if the cell between the new "known"
    cell and the cell required for second order accuracy also happens
    to be "known"), but are not updated in either package.  By luck
    (due to trial values having the same value), the values calculated
    in Scikit-fmm_ for the following example are correct although an
    example that didn't work for Scikit-fmm_ could also be
    constructed.

    >>> phi = computeDistanceFunction([[-1, -1, -1, -1],
    ...                                [ 1,  1, -1, -1],
    ...                                [ 1,  1, -1, -1],
    ...                                [ 1,  1, -1, -1]], order=2)
    >>> phi = computeDistanceFunction(phi, order=2)

    The following values come form Scikit-fmm_.
    
    >>> answer = [[-0.5,        -0.58578644, -1.08578644, -1.85136395],
    ...           [ 0.5,         0.29289322, -0.58578644, -1.54389939],
    ...           [ 1.30473785,  0.5,        -0.5,        -1.5       ],
    ...           [ 1.49547948,  0.5,        -0.5,        -1.5       ]]

    The 3rd and 7th element are different for LSMLIB_. This is because
    the 15th element is not "known" when the "trial" value for the 7th
    element is calculated. Scikit-fmm_ calculates the values in a
    slightly different order so gets a seemingly better answer, but
    this is just chance.

    >>> print np.allclose(phi, answer, rtol=1e-9)
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
    ValueError: dx must be of length len(phi.shape)
    >>> computeExtensionFields([[-1, 1],[1, 1]], speed=[1, 1])
    Traceback (most recent call last):
      ...
    ValueError: phi and speed must have the same shape

    **Test for 1D equality between `distance` and `travel_time`**
    
    >>> phi = np.arange(-5, 5) + 0.499
    >>> d = distance(phi)
    >>> t = travel_time(phi, speed=np.ones_like(phi))
    >>> ##np.testing.assert_allclose(t, np.abs(d))

    **Tests taken from FiPy**

    >>> phi = np.array(((-1, -1, 1, 1),
    ...                 (-1, -1, 1, 1),
    ...                 (1, 1, 1, 1),
    ...                 (1, 1, 1, 1)))
    >>> o1 = distance(phi, order=1)
    >>> dw_o1 =   [[-1.20710678, -0.5,         0.5,         1.5],
    ...            [-0.5,        -0.35355339,  0.5,         1.5],
    ...            [ 0.5,         0.5,         1.20710678,  2.04532893],
    ...            [ 1.5,         1.5,         2.04532893,  2.75243571]]
    >>> np.testing.assert_allclose(o1, dw_o1)

    >>> phi = np.array(((-1, -1, 1, 1),
    ...                 (-1, -1, 1, 1),
    ...                 (1, 1, 1, 1),
    ...                 (1, 1, 1, 1)))
    >>> o1 = travel_time(phi, np.ones_like(phi), order=1)
    >>> dw_o1 =   [[-1.20710678, -0.5,         0.5,         1.5],
    ...            [-0.5,        -0.35355339,  0.5,         1.5],
    ...            [ 0.5,         0.5,         1.20710678,  2.04532893],
    ...            [ 1.5,         1.5,         2.04532893,  2.75243571]]
    >>> ##np.testing.assert_allclose(o1, np.abs(dw_o1))

    >>> phi = np.array(((-1, -1, 1, 1),
    ...                (-1, -1, 1, 1),
    ...                (1, 1, 1, 1),
    ...                (1, 1, 1, 1)))
    >>> o2 = distance(phi)
    >>> dw_o2 = [[-1.30473785,  -0.5,          0.5,         1.49923009],
    ...          [-0.5,        -0.35355339,    0.5,         1.45118446],
    ...          [ 0.5,         0.5,         0.97140452,  1.76215286],
    ...          [ 1.49923009,  1.45118446,  1.76215286,  2.33721352]]

    >>> np.testing.assert_allclose(o2, dw_o2)
    >>> phi = np.array(((-1, -1, 1, 1),
    ...                (-1, -1, 1, 1),
    ...                (1, 1, 1, 1),
    ...                (1, 1, 1, 1)))
    >>> o2 = travel_time(phi, np.ones_like(phi))
    >>> dw_o2 = [[-1.30473785,  -0.5,          0.5,         1.49923009],
    ...          [-0.5,        -0.35355339,    0.5,         1.45118446],
    ...          [ 0.5,         0.5,         0.97140452,  1.76215286],
    ...          [ 1.49923009,  1.45118446,  1.76215286,  2.33721352]]
    >>> ##np.testing.assert_allclose(o2, np.abs(dw_o2))

    >>> distance([-1,1], order=0)
    Traceback (most recent call last):
      ...
    ValueError: order must be 1 or 2

    >>> distance([-1,1], order=3)
    Traceback (most recent call last):
      ...
    ValueError: order must be 1 or 2

    **Extension velocity tests**

    Test 1d extension constant.

    >>> phi =   [-1,-1,-1,1,1,1]
    >>> speed = [1,1,1,1,1,1]
    >>> d, f_ext  = extension_velocities(phi, speed)
    >>> np.testing.assert_allclose(speed, f_ext)

    Test the 1D extension block.

    >>> phi =   np.ones(10)
    >>> phi[0] =- 1
    >>> speed = np.ones(10)
    >>> speed[0:3] = 5
    >>> d, f_ext  = extension_velocities(phi, speed)
    >>> np.testing.assert_allclose(f_ext, 5)

    Test that a uniform speed value is preserved.

    >>> N     = 50
    >>> X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    >>> r     = 0.25
    >>> dx    = 2.0 / (N - 1)
    >>> phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    >>> speed = np.ones_like(phi)
    >>> d, f_ext = extension_velocities(phi, speed, dx)
    >>> np.testing.assert_allclose(f_ext, 1.0)

    Constant value marchout test

    >>> speed[abs(Y)<0.3] = 10.0
    >>> d, f_ext = extension_velocities(phi, speed, dx)
    >>> np.testing.assert_allclose(f_ext, 10.0)

    Test distance from extenstion

    >>> speed = np.ones_like(phi)
    >>> d, f_ext = extension_velocities(phi, speed, dx)
    >>> d2 = distance(phi, dx)
    >>> np.testing.assert_allclose(d, d2)

    Test for extension velocity bug

    >>> N     = 150
    >>> X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    >>> r     = 0.5
    >>> dx    = 2.0 / (N - 1)
    >>> phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    >>> speed = np.ones_like(phi)
    >>> speed[X>0.25] = 3.0
    >>> d2, f_ext = extension_velocities(phi, speed, dx)

    >>> assert (f_ext <= 3.0000001).all()
    >>> assert (f_ext >= 1).all()

    >>> np.testing.assert_almost_equal(f_ext[137, 95], 1, 3)
    >>> np.testing.assert_almost_equal(f_ext[103, 78], 1, 2)
    >>> np.testing.assert_almost_equal(f_ext[72, 100], 3, 3)
    >>> np.testing.assert_almost_equal(f_ext[72, 86], 3, 3)
    >>> np.testing.assert_almost_equal(f_ext[110, 121], 3, 3)

    Simple two point tests

    >>> np.testing.assert_array_equal(distance([-1, 1]),
    ...                                        [-0.5, 0.5])
    >>> np.testing.assert_allclose(distance([-1, -1, -1, 1, 1, 1]),
    ...                                     [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])
    >>> np.testing.assert_allclose(distance([1, 1, 1, -1, -1, -1]),
    ...                                     [2.5, 1.5, 0.5, -0.5, -1.5, -2.5])

    Three point test case

    >>> np.testing.assert_array_almost_equal(distance([-1, 0, 1]),           [-1, 0, 1])
    >>> np.testing.assert_array_almost_equal(distance([-1, 0, 1], dx=[2]),   [-2, 0, 2])
    >>> np.testing.assert_array_almost_equal(distance([-1, 0, 1], dx=2),     [-2, 0, 2])
    >>> np.testing.assert_array_almost_equal(distance([-1, 0, 1], dx=2.0),   [-2, 0, 2])
    >>> np.testing.assert_array_equal(travel_time([1, 0, -1], [1, 1, 1]),
    ...                                           [1, 0, 1])
    >>> np.testing.assert_array_equal(travel_time([-1, 0, 1], [1, 1, 1]),
    ...                                           [1, 0, 1])
    >>> ##np.testing.assert_array_almost_equal(travel_time([1, 0, -1], [1, 1, 1], dx=2), [2, 0, 2])
    >>> ##np.testing.assert_array_equal(travel_time([1, 0, -1], [1, 1, 1], dx=[2]), [2, 0, 2])
    >>> ##np.testing.assert_array_equal(travel_time([1, 0, -1], [1, 1, 1], dx=2.0), [2, 0, 2])

    Travel time tests 1

    >>> ##np.testing.assert_allclose(travel_time([0, 1, 1, 1, 1], [2, 2, 2, 2, 2]), [0, 0.5, 1.0, 1.5, 2.0])
    >>> ##np.testing.assert_array_equal(travel_time([1, 0, -1], [2, 2, 2]), [0.5, 0, 0.5])

    Travel time tests 2
    
    >>> ##phi   = [1, 1, 1, -1, -1, -1]
    >>> ##t     = travel_time(phi, np.ones_like(phi))
    >>> ##exact = [2.5, 1.5, 0.5, 0.5, 1.5, 2.5]
    >>> ##np.testing.assert_allclose(t, exact)

    Travel time tests 3
    
    >>> ##phi   = [-1, -1, -1, 1, 1, 1]
    >>> ##t     = travel_time(phi, np.ones_like(phi))
    >>> ##exact = [2.5, 1.5, 0.5, 0.5, 1.5, 2.5]
    >>> ##np.testing.assert_allclose(t, exact)

    Corner case 

    >>> np.testing.assert_array_almost_equal(distance([0, 0]), [0, 0])
    >>> ##np.testing.assert_array_equal(travel_time([0, 0], [1, 1]), [0, 0])

    Test zero

    >>> distance([1, 0, 1, 1], 0)
    Traceback (most recent call last):
      ...
    ValueError: dx must be greater than zero

    Test dx shape

    >>> distance([0, 0, 1, 0, 0], [0, 0, 1, 0, 0])
    Traceback (most recent call last):
      ...
    ValueError: dx must be of length len(phi.shape)

    Test for small speeds

    Test catching speeds which are too small. Speeds less than the
    machine epsilon are masked off to avoid an overflow.

    >>> ##t = travel_time([-1, -1, 0, 1, 1], [1, 1, 1, 1, 0])
    >>> ##assert isinstance(t, np.ma.MaskedArray)
    >>> ##np.testing.assert_array_equal(t.data[:-1], [2, 1, 0, 1])
    >>> ##np.testing.assert_array_equal(t.mask, [False, False, False, False, True])

    >>> ##t2 = travel_time([-1, -1, 0, 1, 1], [1, 1, 1, 1, 1e-300])
    >>> ##np.testing.assert_array_equal(t, t2)

    Mask test

    Test that when the mask cuts off the solution, the cut off points
    are also masked.

    >>> ##ma    = np.ma.MaskedArray([1, 1, 1, 0], [False, True, False, False])
    >>> ##d     = distance(ma)
    >>> ##exact = np.ma.MaskedArray([0, 0, 1, 0], [True, True, False, False])
    >>> ##np.testing.assert_array_equal(d.mask, exact.mask)
    >>> ##np.testing.assert_array_equal(d, exact)

    Circular level set

    >>> N     = 50
    >>> X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    >>> r     = 0.5
    >>> dx    = 2.0 / (N - 1)
    >>> phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    >>> d     = distance(phi, dx)
    >>> exact = np.sqrt(X ** 2 + Y ** 2) - r
    >>> np.testing.assert_allclose(d, exact, atol=dx)

    Planar level set
    
    >>> N         = 50
    >>> X, Y      = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    >>> dx        = 2.0 / (N - 1)
    >>> phi       = np.ones_like(X)
    >>> phi[0, :] = -1
    >>> d         = distance(phi, dx)
    >>> exact     = Y + 1 - dx / 2.0
    >>> np.testing.assert_allclose(d, exact)

    Masked input
    
    >>> N         = 50
    >>> X, Y      = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    >>> dx        = 2.0 / (N - 1)
    >>> phi       = np.ones_like(X)
    >>> phi[0, 0] = -1
    >>> mask      = np.logical_and(abs(X) < 0.25, abs(Y) < 0.25)
    >>> mphi      = np.ma.MaskedArray(phi.copy(), mask)
    >>> d0        = distance(phi, dx)
    >>> d         = distance(mphi, dx)
    >>> d0[mask]  = 0
    >>> d[mask]   = 0
    >>> shadow    = d0 - d
    >>> bsh       = abs(shadow) > 0.001
    >>> diff      = (bsh).sum()

    >>> ##assert diff > 635 and diff < 645

    Test eikonal solution 

    >>> ##N     = 50
    >>> ##X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    >>> ##r     = 0.5
    >>> ##dx    = 2.0 / (N - 1)
    >>> ##phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    >>> ##speed = np.ones_like(phi) * 2
    >>> ##t     = travel_time(phi, speed, dx)
    >>> ##exact = 0.5 * np.abs(np.sqrt(X ** 2 + Y ** 2) - 0.5)

    >>> ##np.testing.assert_allclose(t, exact, atol=dx)

    Test 1d

    >>> N          = 100
    >>> X          = np.linspace(-1.0, 1.0, N)
    >>> dx         = 2.0 / (N - 1)
    >>> phi        = np.zeros_like(X)
    >>> phi[X < 0] = -1
    >>> phi[X > 0] = 1
    >>> d          = distance(phi, dx)

    >>> np.testing.assert_allclose(d, X)

    Test 3d

    >>> N            = 15
    >>> X            = np.linspace(-1, 1, N)
    >>> Y            = np.linspace(-1, 1, N)
    >>> Z            = np.linspace(-1, 1, N)
    >>> phi          = np.ones((N, N, N))
    >>> phi[0, 0, 0] = -1.0
    >>> dx           = 2.0 / (N - 1)
    >>> d            = distance(phi, dx)
    Traceback (most recent call last):
      ...
    Exception: 3D meshes not yet implemented
    >>> exact        = np.sqrt((X + 1) ** 2 +
    ...                        (Y + 1)[:, np.newaxis] ** 2 +
    ...                        (Z + 1)[:, np.newaxis, np.newaxis] ** 2)

    >>> ## np.testing.assert_allclose(d, exact, atol=dx)

    Test default dx
    
    >>> ##N     = 50
    >>> ##X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    >>> ##r     = 0.5
    >>> ##phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    >>> ##speed = np.ones_like(phi) * 2
    >>> ##out = travel_time(phi, speed)

    Test non-square grid and dx different in different directions
    
    >>> N      = 50
    >>> NX, NY = N, 5 * N
    >>> X, Y   = np.meshgrid(np.linspace(-1, 1, NY), np.linspace(-1, 1, NX))
    >>> r      = 0.5
    >>> phi    = X ** 2 + Y ** 2 - r ** 2
    >>> dx     = [2.0 / (NX - 1), 2.0 / (NY - 1)]
    >>> d      = distance(phi, dx, order=1)
    >>> exact  = np.sqrt(X ** 2 + Y ** 2) - r

    >>> np.testing.assert_allclose(d, exact, atol=1.1*max(dx))
    
    Shape mismatch test

    >>> travel_time([-1, 1], [2])
    Traceback (most recent call last):
      ...
    ValueError: phi and speed must have the same shape

    Speed wrong type test

    >>> travel_time([0, 0, 1, 1], 2)
    Traceback (most recent call last):
      ...
    ValueError: phi and speed must have the same shape

    dx mismatch test

    >>> travel_time([-1, 1], [2, 2], [2, 2, 2, 2])
    Traceback (most recent call last):
      ...
    ValueError: dx must be of length len(phi.shape)

    Test c error handling. Check array type test

    >>> distance(np.array(["a", "b"]))
    Traceback (most recent call last):
      ...
    ValueError: could not convert string to float: a

    
    """

    pass

distance = computeDistanceFunction
extension_velocities = computeExtensionFields

def test():
    r"""
    Run all the doctests available.
    """
    import doctest
    import pylsmlib
    doctest.testmod(pylsmlib)
