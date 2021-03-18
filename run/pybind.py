import pymem3dg as dg

# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/bud.ply"
# refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"
# outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/coarse/Ksg5e-4_H5"
inputMesh = (
    "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//oblate.ply"
)
refMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//patch.ply"
outputDir = (
    "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj"
)

icoFace, icoVertex = dg.getIcosphere(3, 1)
icoFace, icoRefVertex = dg.getIcosphere(3, 1)

# trajFile = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/testrefactor3/traj.nc"

o = dg.Options()
o.isProtein = False
o.isReducedVolume = False
o.isLocalCurvature = True
o.isEdgeFlip = True
o.isGrowMesh = True
o.isVertexShift = False
o.isRefMesh = False
o.isFloatVertex = False
o.isLaplacianMeanCurvature = False

nSub = 1

p = dg.Parameters()
p.Kb = 8.22e-5
p.Kbc = 8.22e-4
p.H0 = -3
p.r_H0 = [0.6, 0.6]
p.eta = 0
p.Ksg = 2
p.Kst = 0
p.Ksl = 1e-7
p.Kse = 1e-7
p.epsilon = -1
p.Bc = -1
p.Kv = 0
p.Vt = -1
p.cam = 0
p.Kf = 0
p.conc = -1
p.height = 0
p.radius = 100000
p.gamma = 0
p.temp = 0
p.pt = [0, 0, 1]
p.lambdaSG = 0
p.lambdaV = 0

h = 1
T = 30000
eps = 0
tSave = 300
isBacktrack = True
rho = 0.99
c1 = 0.0001
isAdaptiveStep = True
verbosity = 3
isAugmentedLagrangian = False
restartNum = 3

p.sigma = (2 * p.gamma * 1.380649e-8 * p.temp / h) ** 0.5

g = dg.System(inputMesh, inputMesh, nSub, p, o)
# g = dg.System(icoFace, icoVertex, icoRefVertex, nSub, p, o)
# g = dg.System(trajFile, -1, nSub, False, p, o)

# g.computeFreeEnergy()
g.computePhysicalForces()

# dg.visualize(g)
# dg.animate_ply(outputDir + "/frames", [1,36], q)

# h = 0.001
# tSave = 0.1
# p.gamma = 0.01
# p.temp = 400
# p.sigma = (2 * p.gamma * 1.380649e-8 * p.temp / h) ** 0.5
# vv = dg.VelocityVerlet(
#     g, h, isAdaptiveStep, T, tSave, eps, outputDir, "/traj.nc", verbosity
# )
# vv.integrate()

# fe = dg.Euler(
#     g,
#     h,
#     isAdaptiveStep,
#     T,
#     tSave,
#     eps,
#     outputDir,
#     "/traj3.nc",
#     verbosity,
#     isBacktrack,
#     rho,
#     c1,
# )
# fe.integrate()

cg = dg.ConjugateGradient(
    g,
    h,
    isAdaptiveStep,
    T,
    tSave,
    eps,
    outputDir,
    "/traj.nc",
    verbosity,
    isBacktrack,
    rho,
    c1,
    0.01,
    isAugmentedLagrangian,
    restartNum,
)
cg.integrate()
# cg.step(1000)
# cg.status()
# cg.saveData()
# cg.march()
# dg.visualize(g)

# bf = dg.BFGS(
#     g,
#     h,
#     isAdaptiveStep,
#     T,
#     tSave,
#     eps,
#     outputDir,
#     "/traj.nc",
#     verbosity,
#     isBacktrack,
#     rho,
#     c1,
#     0.01,
#     isAugmentedLagrangian,
# )
# bf.integrate()

####################################################
# initialize option for visualization, CHANGE HERE #
####################################################
# Quantities data needed for visualization 
Q = dg.Quantities()
Q.ref_coord = True
Q.velocity = True
Q.mean_curvature = True
Q.gauss_curvature = True
Q.spon_curvature = True
Q.ext_pressure = True
Q.physical_pressure = True
Q.capillary_pressure = True
Q.inside_pressure = True
Q.bending_pressure = True
Q.mask = True
Q.H_H0 = True

# GUI & misc: optional arguments for viewers 
# transparency = 1
# angle = 0
# fov = 50
# edgeWidth = 1
# isShow = True
# isSave = False
# screenshotName = "screenshot.png"

# dg.animate_nc(fileName, Q)
