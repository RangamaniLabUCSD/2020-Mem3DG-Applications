import pymem3dg as dg
import numpy as np

####################################################
#                 Initialize pathes                #
####################################################
### Linux ####
inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/results/vesicle/frame100.ply"
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/results/vesicle/frame21.ply"
# refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"
outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/temp"

### Windows ####
# inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//patch.ply"
# inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//frame318.ply"
# refMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//patch.ply"
# refMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//frame24.ply"
# outputDir = (
#     "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//boundaryMutation"
# )

# trajFile = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/testrefactor3/traj.nc"

####################################################
#                 Initialize Meshes                #
####################################################
# patFace, patVertex = dg.getHexagon(1, 4)
# icoFace, icoVertex = dg.getIcosphere(1, 3)
# tetFace, tetVertex = dg.getTetrahedron()
# diaFace, diaVertex = dg.getDiamond(3.14/3)
# cyFace, cyVertex = dg.getCylinder(1, 16, 60, 10.4, 0.1)
# soupFace, soupVertex = dg.processSoup(inputMesh)
####################################################
#        System: Options and Parameters            #
####################################################
o = dg.Options()
o.isReducedVolume = False
o.isConstantOsmoticPressure = True
o.isConstantSurfaceTension = True
o.isProteinVariation = False
o.isShapeVariation = True
o.isFloatVertex = True
o.shapeBoundaryCondition = "fixed"
o.proteinBoundaryCondition = "none"

o.isEdgeFlip = True
o.isSplitEdge = True
o.isCollapseEdge = True
o.isVertexShift = True

nSub = 0
isContinue = False
p = dg.Parameters()
### general ###
p.pt = [0, 0]
p.protein0 = [0.5, 0.5, 1, 0]
### bending ###
p.Kb = 8.22e-5
p.Kbc = 8.22e-5 # 8.22e-4 #DEFINITION OF LARGE AND SMALL VALUE
p.H0c = 6
### surface tension ###
p.Ksg = 1e-3
p.A_res = 0 
p.epsilon = 0
### osmotic force ###
p.Kv = 0
p.V_res = 0
p.Vt = -1
p.cam = -1
### protein binding ###
p.Bc = 0
### line tension ###
p.eta = 0
### DPD ###
p.gamma = 0
p.temp = 0
### regularization ###
p.Kst = 0  # 2e-6
p.Ksl = 0
p.Kse = 0
### Rarely used ###
p.Kf = 0
p.conc = -1
p.height = 0
p.radius = -1
p.lambdaSG = 0
p.lambdaV = 0

# g = dg.System(inputMesh, nSub)
g = dg.System(inputMesh, p, o, nSub, isContinue)
# g = dg.System(soupFace, soupVertex, p, o, nSub)
# g = dg.System(icoFace, icoVertex, p, o, nSub)
# g = dg.System(patFace, patVertex, p, o, nSub)
# g = dg.System(diaFace, diaVertex, diaVertex, nSub, p, o)
# g = dg.System(trajFile, -1, nSub, False, p, o)
# g = dg.System(cyFace, cyVertex, cyVertex, nSub, p, o)
g.meshMutator.flipNonDelaunay = True
g.meshMutator.splitLarge = True
g.meshMutator.splitFat = True
g.meshMutator.splitSkinnyDelaunay = True
g.meshMutator.splitCurved = True
g.meshMutator.curvTol = 0.005
g.meshMutator.collapseSkinny = True
g.computeFreeEnergy()
g.computePhysicalForces()
# dg.visualize(g)

###################################################
#          Time integration / Optimization
####################################################
### options ###
isBacktrack = True
isAdaptiveStep = False
isAugmentedLagrangian = False
### parameters ###
h = 0.05
T = 10000000
eps = 1e-6
tSave = 100
rho = 0.001
c1 = 0.0001
verbosity = 3
restartNum = 3

# h = 0.001
# p.gamma = 0.01
# p.temp = 400
# vv = dg.VelocityVerlet(
#     g, h, isAdaptiveStep, T, tSave, eps, outputDir, "/traj.nc", verbosity
# )
# vv.integrate()

fe = dg.Euler(g, h, T, tSave, eps, outputDir)
fe.tUpdateGeodesics = 50
fe.tProcessMesh = 50
fe.isAdaptiveStep = False
fe.integrate()

# cg = dg.ConjugateGradient(g, h, T, tSave, eps, outputDir)
# cg.tUpdateGeodesics = 50
# cg.tProcessMesh = 50
# cg.restartNum = 1
# cg.isAdaptiveStep = False
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
### Quantities data needed for visualization ###
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