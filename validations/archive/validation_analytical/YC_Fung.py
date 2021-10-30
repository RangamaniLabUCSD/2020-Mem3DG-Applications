import pymem3dg as dg
import polyscope as ps
import numpy as np
from scipy.sparse.linalg import inv

epsilon = 1e-18


def zPrimeDevidedByr(r):
  """dz/dr divided by r

  Args:
      r (numpy array): array of radius

  Returns:
      array of elementwise operation
  """
  return (5614 * r**4 - 10499 * r**2 + 3798) / (2000 * (1-r**2 + epsilon) ** 0.5)


def zPrimePrime(r):
  """dz^2/dr^2

  Args:
      r (numpy array): array of radius

  Returns:
      array of elementwise operation
  """
  return -(22456 * r**6-49068 * r**4+31497 * r**2-3798) / (2000 * (1-r**2 + epsilon)**1.5)


def jacobian_drds(r):
  """jacobian matrix from arclength to radius parametrization, dr/ds

  Args:
      r (numpy array): array of radius

  Returns:
      array of elementwise operation
  """
  return (1 + (zPrimeDevidedByr(r)*r)**2)**(-0.5)

def radius(coord):
  """find X-Y plane radius of the mesh

  Args:
      coord (v times 3 numpy array): mesh vertex matrix

  Returns:
      array of radius
  """
  return (coord[:,0]**2 + coord[:,1]**2)**(0.5)

# c0 = 0.2072
# c1 = 2.0026
# c2 = -1.1228
# R = 1
# coord[:,2] = np.sign(coord[:,2]) * R / 2 * ((1 - r**2)**(0.5)) * (c0 + c1 * r**2 + c2 * r**4)

# deform the icosphere to analytical biconcave shape
topo, coord = dg.getIcosphere(1, 5)
r = radius(coord)
coord[:,2] = np.sign(coord[:,2]) / 2 * ((1 - r**2)**(0.5)) * (0.2072 + 2.0026*r**2 -1.1228*r**4)

# computation of tangential (nu) and trasverse curvature
j = jacobian_drds(r)
kappa_nu = j * zPrimeDevidedByr(r)
kappa_t = j**3 * zPrimePrime(r)
equatorMask = (r==1)
kappa_nu[equatorMask] = -1

# compute analytical and numercial mean curvature 
H = -0.5 * (kappa_t + kappa_nu)
p = dg.Parameters()
p.bending.Kb = 0
p.proteinDistribution.protein0 = [1]
p.bending.Kbc = 1
p.bending.H0c = 0
f = dg.System(topo, coord, p)
f.computeFreeEnergy()
f.computePhysicalForces()
H_num = f.getMeanCurvature()
K_num = f.getGaussianCurvature()
bendingForce = f.forces.getBendingForce()
print(bendingForce)

# compute total mean curvature 
M = f.getLumpedMassMatrix()
totalH_num = np.sum(f.getLumpedMassMatrix() * H_num**2)
totalK_num = np.sum(f.getLumpedMassMatrix() * K_num)
totalH = np.sum(f.getLumpedMassMatrix() * H**2)
volume = f.volume
area = f.surfaceArea
print(volume)
print(area)
print(totalH_num)
print(totalK_num)
print(totalH)
# print(H)

# polyscope visualization
ps.init()
ps_mesh = ps.register_surface_mesh("RBC", coord, topo)
ps_mesh.add_scalar_quantity("mean_curvature1", H, enabled=True)
ps_mesh.add_scalar_quantity("bending force", inv(f.getLumpedMassMatrix()) * np.linalg.norm(bendingForce, axis=1), enabled=True)
ps_mesh.add_scalar_quantity("mean_curvature2", H_num, enabled=True)
ps_mesh.add_scalar_quantity("diff", abs(H + H_num), enabled=True)
ps.set_up_dir("z_up")
ps.show()
