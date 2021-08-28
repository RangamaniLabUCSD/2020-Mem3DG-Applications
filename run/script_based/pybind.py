import pymem3dg as dg
import numpy as np

####################################################
#                 Initialize pathes                #
####################################################
### Linux ####
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/results/temp2/frame1.ply"
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/results/vesicle/frame21.ply"
# refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"
outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/temp2"

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
icoFace, icoVertex = dg.getIcosphere(1, 3)
# tetFace, tetVertex = dg.getTetrahedron()
# diaFace, diaVertex = dg.getDiamond(3.14/3)
# cyFace, cyVertex = dg.getCylinder(1, 16, 60, 7.5, 0)
# soupFace, soupVertex = dg.processSoup(inputMesh)

####################################################
#                 Parameters                       #
####################################################
p = dg.Parameters()
p.proteinMobility = 1
p.temperature = 0
p.point.pt = [0, 0, 1]            
p.point.isFloatVertex = True
p.proteinDistribution.protein0 = [0.1]
p.boundary.shapeBoundaryCondition = "none"
p.boundary.proteinBoundaryCondition = "none"
p.variation.isProteinVariation = True
p.variation.isShapeVariation = True
p.variation.radius = -1
p.bending.Kb = 8.22e-5
p.bending.Kbc = 0  # 8.22e-4 #DEFINITION OF LARGE AND SMALL VALUE
p.bending.H0c = 10
p.tension.isConstantSurfaceTension = False
p.tension.Ksg = 1
p.tension.A_res = 0
p.tension.At = -1
p.tension.lambdaSG = 0
p.adsorption.epsilon = -1e-3
p.osmotic.isPreferredVolume = True
p.osmotic.isConstantOsmoticPressure = False
p.osmotic.Kv = 0.05
p.osmotic.V_res = 0
p.osmotic.n = 1
p.osmotic.Vt = 0.7 * 4.15889  # 1 * 4 * 3.1416 / 3
p.osmotic.cam = -1
p.osmotic.lambdaV = 0
p.dirichlet.eta = 0.1
p.dpd.gamma = 0
p.external.Kf = 0
p.external.conc = -1
p.external.height = 0

####################################################
#                 Mesh processor                   #
####################################################
mP = dg.MeshProcessor()
mP.meshMutator.shiftVertex = False
mP.meshMutator.flipNonDelaunay = True
# mP.meshMutator.splitLarge = True
mP.meshMutator.splitFat = True
mP.meshMutator.splitSkinnyDelaunay = True
# mP.meshMutator.splitCurved = True
# mP.meshMutator.curvTol = 0.005
mP.meshMutator.collapseSkinny = True
# mP.meshRegularizer.Kst = 0.1 # 2e-6
# mP.meshRegularizer.Ksl = 0
# mP.meshRegularizer.Kse = 0
# mP.meshRegularizer.readReferenceData(icoFace, icoVertex, 0)


####################################################
#                 System                           #
####################################################
nSub = 0
isContinue = False
# g = dg.System(inputMesh, nSub)
# g = dg.System(inputMesh, p, nSub, isContinue)
# g = dg.System(soupFace, soupVertex, p, nSub)
g = dg.System(icoFace, icoVertex, p, mP, nSub)
# g = dg.System(patFace, patVertex, p, nSub)
# g = dg.System(diaFace, diaVertex, diaVertex, nSub, p)
# g = dg.System(trajFile, -1, nSub, False, p)
# g = dg.System(cyFace, cyVertex, p, nSub)
# g.computeFreeEnergy()
# g.computePhysicalForces()
# dg.visualize(g)

###################################################
#          Time integration / Optimization
####################################################
### parameters ###
h = 1
T = 10000000
eps = 1e-6
tSave = 300
rho = 0.99
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
# fe.tUpdateGeodesics = 50
fe.tProcessMesh = 300
# fe.isBacktrack = False
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
