import pymem3dg as dg
import numpy as np
from matplotlib import pyplot as plt, ticker as mticker
from scipy.optimize import curve_fit
import polyscope as ps


def myExpFunc(x, a, b):
    return a * np.power(x, b)


faceMatrix, vertexMatrix = dg.getIcosphere(1, 3)
####################################################
#        System: Options and Parameters            #
####################################################
o = dg.Options()
o.isReducedVolume = False
o.isConstantOsmoticPressure = True
o.isConstantSurfaceTension = True
o.isProteinVariation = True
o.isShapeVariation = True
o.isFloatVertex = True
o.shapeBoundaryCondition = "none"
o.proteinBoundaryCondition = "none"

o.isEdgeFlip = False
o.isSplitEdge = False
o.isCollapseEdge = False
o.isVertexShift = False

p = dg.Parameters()
### general ###
p.pt = [0, 0, 1]
p.protein0 = [1, 1, 0.3, 0.2]
p.sharpness = 2
### bending ###
p.Kb = 8.22e-5
p.Kbc = 0
p.H0c = 100
### surface tension ###
p.Ksg = 1e-1
p.A_res = 0
p.epsilon = -1e-1
### osmotic force ###
p.Kv = 1e-1
p.V_res = 0
p.Vt = -1
p.cam = -1
### protein binding ###
p.Bc = 1
### line tension ###
p.eta = 0.1

# ####################################################
# #        System: Options and Parameters            #
# ####################################################
# o = dg.Options()
# o.isReducedVolume = False
# o.isConstantOsmoticPressure = True
# o.isConstantSurfaceTension = True
# o.isProteinVariation = True
# o.isShapeVariation = True
# o.isFloatVertex = True
# o.shapeBoundaryCondition = "none"
# o.proteinBoundaryCondition = "none"

# o.isEdgeFlip = False
# o.isSplitEdge = False
# o.isCollapseEdge = False
# o.isVertexShift = False

# p = dg.Parameters()
# ### general ###
# p.pt = [0, 0, 1]
# p.protein0 = [1, 1, 0.3, 0.1]
# ### bending ###
# p.Kb = 8.22e-5
# p.Kbc = 0  # 8.22e-4 #DEFINITION OF LARGE AND SMALL VALUE
# p.H0c = 10
# ### surface tension ###
# p.Ksg = 1e-3
# p.A_res = 0
# p.epsilon = -1e-3
# ### osmotic force ###
# p.Kv = 1e-3
# p.V_res = 0
# p.Vt = -1
# p.cam = -1
# ### protein binding ###
# p.Bc = 1
# ### line tension ###
# p.eta = 0.01


nSubSize = 8
maxh = 0.1

H = np.zeros(nSubSize)
diffBendingForce = np.zeros(nSubSize)
diffCapillaryForce = np.zeros(nSubSize)
diffOsmoticForce = np.zeros(nSubSize)
diffLineCapillaryForce = np.zeros(nSubSize)
diffAdsorptionForce = np.zeros(nSubSize)
diffBendingPotential = np.zeros(nSubSize)
diffAdsorptionPotential = np.zeros(nSubSize)
diffDiffusionPotential = np.zeros(nSubSize)

for i in range(nSubSize):
    h = maxh / (2**i)
    H[i] = h

    g = dg.System(faceMatrix, vertexMatrix, p, o)
    g.computeFreeEnergy()
    g.computePhysicalForces()
    p.protein0 = g.getProteinDensity()
    # print("new pos: ", np.linalg.norm(h * g.F.getBendingForce() +
    #                  g.getVertexPositionMatrix()))

    gNew = dg.System(faceMatrix, h * g.F.getBendingForce() +
                     vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.BE - gNew.E.BE
    expected = h * np.linalg.norm(g.F.getBendingForce())**2
    diffBendingForce[i] = abs(dEnergy - expected)
    print("h is", h, " and dBE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.F.getCapillaryForce() +
                     vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.sE - gNew.E.sE
    expected = h * np.linalg.norm(g.F.getCapillaryForce())**2
    diffCapillaryForce[i] = abs(dEnergy - expected)
    print("h is", h, " and dsE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.F.getOsmoticForce() +
                     vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.pE - gNew.E.pE
    expected = h * np.linalg.norm(g.F.getOsmoticForce())**2
    diffOsmoticForce[i] = abs(dEnergy - expected)
    print("h is", h, " and dpE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.F.getAdsorptionForce() +
                     vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.aE - gNew.E.aE
    expected = h * np.linalg.norm(g.F.getAdsorptionForce())**2
    diffAdsorptionForce[i] = abs(dEnergy - expected)
    print("h is", h, " and daE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.F.getLineCapillaryForce() +
                     vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.dE - gNew.E.dE
    expected = h * np.linalg.norm(g.F.getLineCapillaryForce())**2
    diffLineCapillaryForce[i] = abs(dEnergy - expected)
    print("h is", h, " and ddE: ", dEnergy, " and expected: ", expected)

    p.protein0 = g.getProteinDensity() + h * g.F.getDiffusionPotential()
    gNew = dg.System(faceMatrix, vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.dE - gNew.E.dE
    expected = h * np.linalg.norm(g.F.getDiffusionPotential())**2
    diffDiffusionPotential[i] = abs(dEnergy - expected)
    print("h is", h, " and ddE: ", dEnergy, " and expected: ", expected)

    p.protein0 = g.getProteinDensity() + h * g.F.getBendingPotential()
    gNew = dg.System(faceMatrix, vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.BE - gNew.E.BE
    expected = h * np.linalg.norm(g.F.getBendingPotential())**2
    diffBendingPotential[i] = abs(dEnergy - expected)
    print("h is", h, " and dBE: ", dEnergy, " and expected: ", expected)

    p.protein0 = g.getProteinDensity() + h * g.F.getAdsorptionPotential()
    gNew = dg.System(faceMatrix, vertexMatrix, p, o)
    gNew.computeFreeEnergy()
    dEnergy = g.E.aE - gNew.E.aE
    expected = h * np.linalg.norm(g.F.getAdsorptionPotential())**2
    diffAdsorptionPotential[i] = abs(dEnergy - expected)
    print("h is", h, " and daE: ", dEnergy, " and expected: ", expected)


plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
SMALL_SIZE = 7
MEDIUM_SIZE = 9
BIGGER_SIZE = 10
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('pdf', fonttype=42)
# plt.subplots_adjust(left=0.164, bottom=0.07, right=0.988, top=0.988)

fig, axs = plt.subplots(1, 3)
fig.set_size_inches(7, 4)

popt, pcov = curve_fit(myExpFunc, H, diffBendingForce)
axs[0].loglog(H, diffBendingForce,
              '-o',  label="$f^b, O({{{order: .1f}}})$".format(order=popt[1]))

popt, pcov = curve_fit(myExpFunc, H, diffCapillaryForce)
axs[0].loglog(H, diffCapillaryForce,
              '-o',  label="$f^s, O({{{order: .1f}}})$".format(order=popt[1]))


popt, pcov = curve_fit(myExpFunc, H, diffOsmoticForce)
axs[0].loglog(H, diffOsmoticForce,
              '-o',  label="$f^p, O({{{order: .1f}}})$".format(order=popt[1]))


popt, pcov = curve_fit(myExpFunc, H, diffLineCapillaryForce)
axs[0].loglog(H, diffLineCapillaryForce,
              '-o',  label="$f^d, O({{{order: .1f}}})$".format(order=popt[1]))

popt, pcov = curve_fit(myExpFunc, H, diffAdsorptionForce)
axs[0].loglog(H, diffAdsorptionForce,
              '-o',  label="$f^a, O({{{order: .1f}}})$".format(order=popt[1]))

popt, pcov = curve_fit(myExpFunc, H, diffDiffusionPotential)
axs[0].loglog(H, diffDiffusionPotential,
              '-o',  label="$\mu^d, O({{{order: .1f}}})$".format(order=popt[1]))

popt, pcov = curve_fit(myExpFunc, H, diffBendingPotential)
axs[0].loglog(H, diffBendingPotential,
              '-o',  label="$\mu^b, O({{{order: .1f}}})$".format(order=popt[1]))

print("adsorption error: ", diffAdsorptionPotential)
popt, pcov = curve_fit(myExpFunc, H, diffAdsorptionPotential)
axs[0].loglog(H, diffAdsorptionPotential,
              '-o',  label="$\mu^a, O({{{order: .1f}}})$".format(order=popt[1]))

axs[0].legend()
plt.tight_layout()
plt.savefig("gradient.pdf", transparent=True)

plt.show()
