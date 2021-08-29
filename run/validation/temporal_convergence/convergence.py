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
p = dg.Parameters()
p.point.isFloatVertex = True
p.point.pt = [0, 0, 1]
p.boundary.shapeBoundaryCondition = "none"
p.boundary.proteinBoundaryCondition = "none"
p.variation.isProteinVariation = True
p.variation.isShapeVariation = True
p.proteinDistribution.protein0 = [1, 1, 0.3, 0.2]
p.proteinDistribution.sharpness = 2
p.bending.Kb = 8.22e-5
p.bending.Kbc = 0
p.bending.H0c = 1
p.tension.isConstantSurfaceTension = True
p.tension.Ksg = 1e-1
p.adsorption.epsilon = -1e-1
p.osmotic.isConstantOsmoticPressure = True
p.osmotic.Kv = 1e-1
p.proteinMobility = 1
p.dirichlet.eta = 0.1


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

    g = dg.System(faceMatrix, vertexMatrix, p)
    g.computeFreeEnergy()
    g.computePhysicalForces()
    p.proteinDistribution.protein0 = g.getProteinDensity()
    # print("new pos: ", np.linalg.norm(h * g.forces.getBendingForce() +
    #                  g.getVertexPositionMatrix()))

    gNew = dg.System(faceMatrix, h * g.forces.getBendingForce() +
                     vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.BE - gNew.energy.BE
    expected = h * np.linalg.norm(g.forces.getBendingForce())**2
    diffBendingForce[i] = abs(dEnergy - expected)
    print("h is", h, " and dBE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.forces.getCapillaryForce() +
                     vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.sE - gNew.energy.sE
    expected = h * np.linalg.norm(g.forces.getCapillaryForce())**2
    diffCapillaryForce[i] = abs(dEnergy - expected)
    print("h is", h, " and dsE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.forces.getOsmoticForce() +
                     vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.pE - gNew.energy.pE
    expected = h * np.linalg.norm(g.forces.getOsmoticForce())**2
    diffOsmoticForce[i] = abs(dEnergy - expected)
    print("h is", h, " and dpE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.forces.getAdsorptionForce() +
                     vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.aE - gNew.energy.aE
    expected = h * np.linalg.norm(g.forces.getAdsorptionForce())**2
    diffAdsorptionForce[i] = abs(dEnergy - expected)
    print("h is", h, " and daE: ", dEnergy, " and expected: ", expected)

    gNew = dg.System(faceMatrix, h * g.forces.getLineCapillaryForce() +
                     vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.dE - gNew.energy.dE
    expected = h * np.linalg.norm(g.forces.getLineCapillaryForce())**2
    diffLineCapillaryForce[i] = abs(dEnergy - expected)
    print("h is", h, " and ddE: ", dEnergy, " and expected: ", expected)

    p.proteinDistribution.protein0 = g.getProteinDensity() + h * g.forces.getDiffusionPotential()
    gNew = dg.System(faceMatrix, vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.dE - gNew.energy.dE
    expected = h * np.linalg.norm(g.forces.getDiffusionPotential())**2
    diffDiffusionPotential[i] = abs(dEnergy - expected)
    print("h is", h, " and ddE: ", dEnergy, " and expected: ", expected)

    p.proteinDistribution.protein0 = g.getProteinDensity() + h * g.forces.getBendingPotential()
    gNew = dg.System(faceMatrix, vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.BE - gNew.energy.BE
    expected = h * np.linalg.norm(g.forces.getBendingPotential())**2
    diffBendingPotential[i] = abs(dEnergy - expected)
    print("h is", h, " and dBE: ", dEnergy, " and expected: ", expected)

    p.proteinDistribution.protein0 = g.getProteinDensity() + h * g.forces.getAdsorptionPotential()
    gNew = dg.System(faceMatrix, vertexMatrix, p)
    gNew.computeFreeEnergy()
    dEnergy = g.energy.aE - gNew.energy.aE
    expected = h * np.linalg.norm(g.forces.getAdsorptionPotential())**2
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
