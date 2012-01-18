from fipy import CellVariable, Viewer, Grid2D
import pyximport
import numpy as np
pyximport.install(setup_args = {'options' :
                                {'build_ext' :
                                 {'libraries' : ['lsm_serial', 'lsm_toolbox'],
                                  'include_dirs' : np.get_include()
                                  }}})


from pylsmlib import computeDistanceFunction2d_

N = 100
L = 1.

m = Grid2D(nx=N, ny=N, Lx=L, Ly=L)
phi = CellVariable(mesh=m)
phi[:] = (m.x - L / 2)**2 + (m.y - L / 2)**2 - 0.25**2
vi = Viewer(phi)
vi.plot()
phi[:] = computeDistanceFunction2d_(phi.value, nx=m.nx, ny=m.ny, dx=m.dx, dy=m.dy)
vi.plot()
raw_input('stopped')

