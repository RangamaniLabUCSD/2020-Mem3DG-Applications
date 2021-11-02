import pymem3dg as dg
import numpy as np
import parameters

####################################################
#                 Initialize pathes                #
####################################################
""" Linux """
outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/temp3"

""" Windows """
# outputDir = (
#     "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//bud//asymm//testTraj//boundaryMutation"
# )

# trajFile = "/home/cuzhu/2020-Mem3DG-Applications/results/bud/testrefactor3/traj.nc"

####################################################
#            Initialize input geometry             #
####################################################
""" Built-in construction """
# patFace, patVertex = dg.getHexagon(1, 4)
icoFace, icoVertex = dg.getIcosphere(1, 3)
# tetFace, tetVertex = dg.getTetrahedron()
# diaFace, diaVertex = dg.getDiamond(3.14/3)
# cyFace, cyVertex = dg.getCylinder(1, 16, 60, 7.5, 0)
# soupFace, soupVertex = dg.processSoup(inputMesh)

""" Linux """
# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/patch.ply"

""" Windows """
# inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//patch.ply"

####################################################
#                 Parameters                       #
####################################################
""" Import from file """
p = parameters.parameters()

####################################################
#                 Mesh processor                   #
####################################################
mP = dg.MeshProcessor()
mP.meshMutator.shiftVertex = True
mP.meshMutator.flipNonDelaunay = True
# mP.meshMutator.splitLarge = True
mP.meshMutator.splitFat = True
mP.meshMutator.splitSkinnyDelaunay = True
mP.meshMutator.splitCurved = True
mP.meshMutator.curvTol = 0.005
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

""" System construction """
# g = dg.System(inputMesh, nSub)
# g = dg.System(inputMesh, p, nSub, isContinue)
# g = dg.System(soupFace, soupVertex, p, nSub)
g = dg.System(icoFace, icoVertex, p, mP, nSub)
# g = dg.System(patFace, patVertex, p, nSub)
# g = dg.System(diaFace, diaVertex, diaVertex, nSub, p)
# g = dg.System(trajFile, -1, nSub, False, p)
# g = dg.System(cyFace, cyVertex, p, nSub)

###################################################
#          Time integration / Optimization
####################################################
""" Integrator setups (essential) """
h = 1
T = 10000000
eps = 1e-6
tSave = 300
verbosity = 5

""" Integrator construction """
fe = dg.Euler(g, h, T, tSave, eps, outputDir)

""" Integrator setups (optional) """
# fe.tUpdateGeodesics = 50
fe.processMeshPeriod = 50
# fe.isBacktrack = False
fe.isAdaptiveStep = True
fe.verbosity = verbosity
fe.integrate()