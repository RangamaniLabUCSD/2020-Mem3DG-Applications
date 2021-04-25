import pymem3dg as dg

####################################################
#                 Initialize pathes                #
####################################################
### Linux ####
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/results/bud//asymm/testTraj/frame320.ply"
# refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"
# outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/asymm/testTraj"

### Windows ####
# inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//patch.ply"
inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//frame318.ply"
# refMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//patch.ply"
# refMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//frame24.ply"
outputDir = (
    "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//boundaryMutation"
)

# trajFile = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/testrefactor3/traj.nc"

####################################################
#                 Initialize Meshes                #
####################################################
# icoFace, icoVertex = dg.getIcosphere(0, 1)
# icoFace, icoRefVertex = dg.getIcosphere(0, 1)
# tetFace, tetVertex = dg.getTetrahedron()
# diaFace, diaVertex = dg.getDiamond(3.14/3)


####################################################
#        System: Options and Parameters            #
####################################################
o = dg.Options()
o.isProtein = False
o.isReducedVolume = False
o.isLocalCurvature = True
o.isEdgeFlip = True
o.isGrowMesh = True
o.isVertexShift = False
o.isRefMesh = False
o.isFloatVertex = True
o.isLaplacianMeanCurvature = False

nSub = 0

p = dg.Parameters()
p.Kb = 8.22e-5
p.Kbc = 8.22e-4  # 8.22e-4 DEFINITION OF LARGE AND SMALL VALUE
p.H0 = 6
p.r_H0 = [0.5, 0.5]
p.eta = 0
p.Ksg = 1e-2  # 1e-2 DEFINITION OF LARGE AND SMALL VALUE
p.Kst = 0  # 2e-6
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
p.pt = [0, 0]
p.lambdaSG = 0
p.lambdaV = 0

g = dg.System(inputMesh, inputMesh, nSub, p, o)
# g = dg.System(icoFace, icoVertex, icoRefVertex, nSub, p, o)
# g = dg.System(diaFace, diaVertex, diaVertex, nSub, p, o)
# g = dg.System(trajFile, -1, nSub, False, p, o)

# g.computeFreeEnergy()
# g.computePhysicalForces()
# dg.visualize(g)

####################################################
#          Time integration / Optimization         #
####################################################
h = 0.1
T = 5000
eps = 0
tSave = 10
isBacktrack = True
rho = 0.99
c1 = 0.0001
isAdaptiveStep = True
verbosity = 3
isAugmentedLagrangian = False
restartNum = 3

# h = 0.001
# p.gamma = 0.01
# p.temp = 400
# vv = dg.VelocityVerlet(
#     g, h, isAdaptiveStep, T, tSave, eps, outputDir, "/traj.nc", verbosity
# )
# vv.integrate()

fe = dg.Euler(
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
)
fe.integrate()

# cg = dg.ConjugateGradient(
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
#     restartNum,
# )
# cg.integrate()
# cg.step(1000)
# cg.status()
# cg.saveData()
# cg.march()

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
Q.ext_force = True
Q.physical_force = True
Q.capillary_force = True
Q.osmotic_force = True
Q.bending_force = True
Q.mask = True
Q.H_H0 = True
Q.the_point = True

# GUI & misc: optional arguments for viewers
# transparency = 1
# angle = 0
# fov = 50
# edgeWidth = 1
# isShow = True
# isSave = False
# screenshotName = "screenshot.png"

# dg.animate_nc(fileName, Q)
# dg.animate_ply(outputDir + "/frames", [1,36], q)
