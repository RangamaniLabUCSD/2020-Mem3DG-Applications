import pymem3dg

inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/oblate.ply"
refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/oblate.ply"
outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/testrefactor4"
# inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//slightlyOblate.ply"
# refMesh =  "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//slightlyOblate.ply"
# outputDir = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//vesicle//testAPI"

trajFile = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/testrefactor3/traj.nc"

isProtein = False
isVertexShift = False
isReducedVolume = False
isLocalCurvature = False

nSub = 0

p = pymem3dg.Parameters()
p.Kb = 8.22e-5
p.H0 = 0
p.r_H0 = [-1, -1]
p.eta = 0
p.Ksg = 1
p.Kst = 1e-8 #2e-6    
p.Ksl = 0
p.Kse = 0
p.epsilon = -1
p.Bc = -1
p.Kv = 0.1
p.Vt = -1
p.cam = 0.32
p.Kf = 0
p.conc = -1
p.height = 0
p.radius = 100000
p.gamma = 0
p.temp = 0
p.pt = [0, 0, 0]
p.lambdaSG = 0
p.lambdaV = 0

h = 0.1
T = 100000
eps = 0.001
tSave = 1
isBacktrack =  True
rho = 0.99
c1 = 0.0001
isAdaptiveStep = True
verbosity = 3
isAugmentedLagrangian = False

p.sigma = (2 * p.gamma * 1.380649e-8 * p.temp / h)**0.5


isLocalCurvature=True
p.r_H0=[1, 1]
p.eta=0.0005
p.H0 = 3
g = pymem3dg.System(inputMesh, refMesh, nSub, p, isReducedVolume,
                    isProtein, isLocalCurvature, isVertexShift)
# g = pymem3dg.System(trajFile, -1, nSub, False, p,
#                     isReducedVolume, isProtein, isLocalCurvature, isVertexShift)

# g.computeFreeEnergy()
# g.computeAllForces()
# g.visualize()
# fb = f.getBendingPressure()
# totalEnergy = f.E.totalE

# h = 0.001
# tSave = 0.1
# p.gamma = 0.01;
# p.temp = 400
# p.sigma = (2 * p.gamma * 1.380649e-8 * p.temp / h)**0.5
# bf = pymem3dg.VelocityVerlet(g, h, isAdaptiveStep, T,
#                     tSave, eps, outputDir, "/traj.nc", verbosity)
# bf.integrate()

# fe = pymem3dg.Euler(g, h, isAdaptiveStep, T,
#                     tSave, eps, outputDir, "/traj3.nc", verbosity, isBacktrack, rho, c1)
# fe.integrate()

bf = pymem3dg.ConjugateGradient(g, h, isAdaptiveStep, T,
                    tSave, eps, outputDir, "/traj.nc", verbosity, isBacktrack, rho, c1, 0.01, isAugmentedLagrangian)
bf.integrate()

# cg.step(1000)
# cg.status()
# cg.saveData()
# cg.march()

# pymem3dg.animation_nc(fileName=outputDir + "/traj3.nc", transparency=0.5, angle=0,
#                       fov=50, edgeWidth=2, ref_coord=True, velocity=True,
#                       mean_curvature=True, gauss_curvature=True, spon_curvature=True,
#                       ext_pressure=True, physical_pressure=True,
#                       capillary_pressure=True, inside_pressure=True,
#                       bending_pressure=True, line_pressure=True, mask=True, H_H0=True)
