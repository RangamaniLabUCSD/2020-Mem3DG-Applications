import pymem3dg

inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/slightlyOblate.ply"
refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/slightlyOblate.ply"
outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/vesicle/testAPI"
# inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//slightlyOblate.ply"
# refMesh =  "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//slightlyOblate.ply"
# outputDir = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//vesicle//testAPI"

isProtein = False
isVertexShift = False
isReducedVolume = False
isLocalCurvature = True

nSub = 0

p = pymem3dg.Parameters()
p.Kb = 8.22e-5
p.H0 = 3
p.sharpness = 50
p.r_H0 = [1, 1]
p.eta = 0.002
p.Ksg = 1
p.Kst = 0.1
p.Ksl = 0
p.Kse = 0
p.epsilon = -1
p.Bc = -1
p.Kv = 0
p.Vt = -1
p.cam = 0.3
p.Kf = 0
p.conc = -1
p.height = 0
p.radius = 100000
p.gamma = 0
p.temp = 0
p.pt = [0, 0, 0]
p.lambdaSG = 0
p.lambdaV = 0

h = 0.005
T = 1000
eps = 1e-6
tSave = 10
isBacktrack = True
rho = 0.5
c1 = 0
isAdaptiveStep = True
verbosity = 3
isAugmentedLagrangian = False

p.sigma = (2 * p.gamma * 1.380649e-8 * p.temp / h)**0.5

g = pymem3dg.System(inputMesh, refMesh, nSub, p, isReducedVolume,
                    isProtein, isLocalCurvature, isVertexShift)

g.computeFreeEnergy()
print(g.E.totalE)
print(g.E.lE)
g.computeAllForces()
print(g.getLineTensionPressure())
# g.visualize()

fb = f.getBendingPressure()
totalEnergy = f.E.totalE

vv = pymem3dg.VelocityVerlet(f, h, isAdaptiveStep, T,
                             tSave, eps, outputDir, "/traj3.nc", verbosity)
vv.integrate()

fe = pymem3dg.Euler(g, h, isAdaptiveStep, T,
                    tSave, eps, outputDir, "/traj3.nc", verbosity, isBacktrack, rho, c1)
fe.integrate();

pymem3dg.animation_nc(fileName=outputDir + "/traj3.nc", transparency=0.5, angle=0,
                              fov=50, edgeWidth=2, ref_coord=True, velocity=True,
                              mean_curvature=True, gauss_curvature=True, spon_curvature=True,
                              ext_pressure=True, physical_pressure=True,
                              capillary_pressure=True, inside_pressure=True,
                              bending_pressure=True, line_pressure=True, mask=True, H_H0=True)
fe.step(1000)
g.visualize()

cg = pymem3dg.ConjugateGradient(f, h, isAdaptiveStep, T,
                                tSave, eps, outputDir, "/traj.nc", verbosity, isBacktrack, rho, c1, 0.01, isAugmentedLagrangian)
cg.integrate()
print(f.E.potE)
cg.status()
cg.saveData()
cg.march()
print(f.E.potE)
cg.step(100)
print(f.E.potE)
