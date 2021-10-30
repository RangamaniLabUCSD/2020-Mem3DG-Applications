import pymem3dg as dg
import numpy as np
from matplotlib import pyplot as plt, ticker as mticker
from scipy.optimize import curve_fit
import polyscope as ps
####################################################
#        System: Options and Parameters            #
####################################################
p = dg.Parameters()
p.point.isFloatVertex = True
p.point.pt = [0, 0, 1]
p.boundary.shapeBoundaryCondition = "none"
p.boundary.proteinBoundaryCondition = "none"
p.variation.isProteinVariation = False
p.variation.isShapeVariation = True
p.proteinDistribution.protein0 = [1]
p.proteinDistribution.sharpness = 2
p.bending.Kb = 0
p.bending.Kbc = 8.22e-5
p.bending.H0c = 1.5
p.tension.isConstantSurfaceTension = True
p.tension.Ksg = 1e-1
p.adsorption.epsilon = 0
p.osmotic.isConstantOsmoticPressure = True
p.osmotic.Kv = 1e-1
p.proteinMobility = 0
p.dirichlet.eta = 0


def myExpFunc(x, a, b):
    return a * np.power(x, b)


# Set unit of data
kBT = 1.38e-23 * 310 * 1e15  # 1e-15 J
pJ = 1e3  # fJ = 1e-15 J
Pa = 1e-3  # kPa

nSubSize = 7
minSub = 1
nsubdivisions = np.zeros(nSubSize)
nvertices = np.zeros(nSubSize)
normBendingForce = np.zeros(nSubSize)
normBendingForce_areaGrad = np.zeros(nSubSize)
normBendingForce_gaussVec = np.zeros(nSubSize)
normBendingForce_schlafliVec = np.zeros(nSubSize)
normCapillaryForce = np.zeros(nSubSize)
normOsmoticForce = np.zeros(nSubSize)
# normLineCapillaryForce = np.zeros(nSubSize)
# normAdsorptionForce = np.zeros(nSubSize)
# normBendingPotential = np.zeros(nSubSize)
# normAdsorptionPotential = np.zeros(nSubSize)
# normDiffusionPotential = np.zeros(nSubSize)

BE = np.zeros(nSubSize)
sE = np.zeros(nSubSize)
pE = np.zeros(nSubSize)
# dE = np.zeros(nSubSize)
# aE = np.zeros(nSubSize)

for i in range(nSubSize):
    n = minSub + i
    faceMatrix, vertexMatrix = dg.getIcosphere(1, n)
    g = dg.System(faceMatrix, vertexMatrix, p)
    g.computeTotalEnergy()
    g.computePhysicalForces()
    g.saveRichData("nSub" + str(n) + ".ply")

    nvertices[i] = len(np.array(g.getVertexPositionMatrix()))
    nsubdivisions[i] = n

    # L1 norm
    normBendingForce[i] = np.sum(abs(g.forces.getBendingForce()))
    normBendingForce_areaGrad[i] = np.sum(abs(g.forces.getBendingForce_areaGrad()))
    normBendingForce_gaussVec[i] = np.sum(abs(g.forces.getBendingForce_gaussVec()))
    normBendingForce_schlafliVec[i] = np.sum(
        abs(g.forces.getBendingForce_schlafliVec()))
    normCapillaryForce[i] = np.sum(abs(g.forces.getCapillaryForce()))
    normOsmoticForce[i] = np.sum(abs(g.forces.getOsmoticForce()))
    # # normLineCapillaryForce[i] = np.sum(abs(g.forces.getLineCapillaryForce()))
    # normAdsorptionForce[i] = np.sum(abs(g.forces.getAdsorptionForce()))
    # normBendingPotential[i] = np.sum(abs(g.forces.getBendingPotential()))
    # normDiffusionPotential[i] = np.sum(abs(g.forces.getDiffusionPotential()))
    # normAdsorptionPotential[i] = np.sum(abs(g.forces.getAdsorptionPotential()))

    # summation
    # normBendingForce[i] = np.sum(g.forces.getBendingForce())
    # normCapillaryForce[i] = np.sum(g.forces.getCapillaryForce())
    # normOsmoticForce[i] = np.sum(g.forces.getOsmoticForce())
    # normLineCapillaryForce[i] = np.sum(g.forces.getLineCapillaryForce())
    # normAdsorptionForce[i] = np.sum(g.forces.getAdsorptionForce())

    # L2 norm
    # normBendingForce[i] = np.linalg.norm(g.forces.getBendingForce())
    # normCapillaryForce[i] = np.linalg.norm(g.forces.getCapillaryForce())
    # normOsmoticForce[i] = np.linalg.norm(g.forces.getOsmoticForce())
    # normLineCapillaryForce[i] = np.linalg.norm(g.forces.getLineCapillaryForce())
    # normAdsorptionForce[i] = np.linalg.norm(g.forces.getAdsorptionForce())

    BE[i] = g.energy.bendingEnergy
    sE[i] = g.energy.surfaceEnergy
    pE[i] = g.energy.pressureEnergy
    # dE[i] = g.energy.dE
    # aE[i] = g.energy.aE

print("normBendingForce: ", normBendingForce)
# # print("normLineCapillaryForce: ", normLineCapillaryForce)
print("normOsmoticForce: ", normOsmoticForce)
print("normCapillaryForce: ", normCapillaryForce)
print("BE: ", BE)
print("sE: ", sE)
print("pE: ", pE)
# print("lE: ", dE)
# Visualization preference
# Font:
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
# mpl.rcParams.update({'font.size': 8})
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

# x = np.power(nvertices, -0.5)[:nSubSize-1]
# x = 1 / 2 ** np.array(range(nSubSize-1))
x = 2 ** np.array(list(reversed(range(nSubSize-1))))

error = abs(normCapillaryForce - normCapillaryForce[nSubSize-1]
            )[:nSubSize-1] / abs(normCapillaryForce[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[0].loglog(x, error,
              '-o', label="$f^s, O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(normOsmoticForce - normOsmoticForce[nSubSize-1]
            )[:nSubSize-1] / abs(normOsmoticForce[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[0].loglog(x, error,
              '-o', label="$f^p, O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(normBendingForce - normBendingForce[nSubSize-1]
            )[:nSubSize-1] / abs(normBendingForce[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[0].loglog(x, error,
              '-o',  label="$f^b, O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(normBendingForce_areaGrad - normBendingForce_areaGrad[nSubSize-1]
            )[:nSubSize-1] / abs(normBendingForce_areaGrad[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[0].loglog(x, error,
              '--',  label="$f^b (H), O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(normBendingForce_gaussVec - normBendingForce_gaussVec[nSubSize-1]
            )[:nSubSize-1] / abs(normBendingForce_gaussVec[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[0].loglog(x, error,
              '--',  label="$f^b (K), O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(normBendingForce_schlafliVec - normBendingForce_schlafliVec[nSubSize-1]
            )[:nSubSize-1] / abs(normBendingForce_schlafliVec[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[0].loglog(x, error,
              '--',  label="$f^b (S), O({{{order: .1f}}})$".format(order=popt[1]))

# # error = abs(normLineCapillaryForce - normLineCapillaryForce[nSubSize-1]
# )[:nSubSize-1] / abs(normLineCapillaryForce[nSubSize-1])
# popt, pcov = curve_fit(myExpFunc, x, error)
# axs[0].loglog(x, error,
#               '-o', label="$f^d, O({{{order: .1f}}})$".format(order=popt[1]))

# error = abs(normAdsorptionForce - normAdsorptionForce[nSubSize-1]
#             )[:nSubSize-1] / abs(normAdsorptionForce[nSubSize-1])
# popt, pcov = curve_fit(myExpFunc, x, error)
# axs[0].loglog(x, error,
#               '-o', label="$f^a, O({{{order: .1f}}})$".format(order=popt[1]))

# error = abs(normAdsorptionPotential - normAdsorptionPotential[nSubSize-1]
#             )[:nSubSize-1] / abs(normAdsorptionPotential[nSubSize-1])
# popt, pcov = curve_fit(myExpFunc, x, error)
# axs[0].loglog(x, error,
#               '-o', label="$\mu^a, O({{{order: .1f}}})$".format(order=popt[1]))

# # error = abs(normDiffusionPotential - normDiffusionPotential[nSubSize-1]
# )[:nSubSize-1] / abs(normDiffusionPotential[nSubSize-1])
# popt, pcov = curve_fit(myExpFunc, x, error)
# axs[0].loglog(x, error,
#               '-o', label="$\mu^d, O({{{order: .1f}}})$".format(order=popt[1]))

# error = abs(normBendingPotential - normBendingPotential[nSubSize-1]
#             )[:nSubSize-1] / abs(normBendingPotential[nSubSize-1])
# popt, pcov = curve_fit(myExpFunc, x, error)
# axs[0].loglog(x, error,
#               '-o', label="$\mu^b, O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(sE - sE[nSubSize-1]
            )[:nSubSize-1] / abs(sE[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[1].loglog(x, error,
              '-o', label="$E_s, O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(pE - pE[nSubSize-1]
            )[:nSubSize-1] / abs(pE[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[1].loglog(x, error,
              '-o', label="$E_p, O({{{order: .1f}}})$".format(order=popt[1]))

error = abs(BE - BE[nSubSize-1]
            )[:nSubSize-1] / abs(BE[nSubSize-1])
popt, pcov = curve_fit(myExpFunc, x, error)
axs[1].loglog(x, error,
              '-o', label="$E_b, O({{{order: .1f}}})$".format(order=popt[1]))

# error = abs(aE - aE[nSubSize-1]
#             )[:nSubSize-1] / abs(aE[nSubSize-1])
# popt, pcov = curve_fit(myExpFunc, x, error)
# axs[1].loglog(x, error,
#               '-o', label="$E_a, O({{{order: .1f}}})$".format(order=popt[1]))

# # error = abs(dE - dE[nSubSize-1]
# )[:nSubSize-1] / abs(dE[nSubSize-1])
# popt, pcov = curve_fit(myExpFunc, x, error)
# axs[1].loglog(x, error,
#               '-o', label="$E_d, O({{{order: .1f}}})$".format(order=popt[1]))

axs[0].legend()
axs[1].legend()
# axs[0].xaxis.set_minor_formatter(mticker.ScalarFormatter())
axs[0].set_xticks([1, 10])
axs[1].set_xticks([1, 10])


# ax.legend()
# ax.set_ylabel("$e_E / E$")
# ax.set_xlabel("$L$")
plt.tight_layout()
plt.savefig("homo.pdf", transparent=True)

plt.show()


# outputMesh = "nSub2.ply"

# ps.init()
# face, vertex = dg.readMesh(outputMesh)
# mesh = ps.register_surface_mesh("frame0", vertex, face)
# mesh.add_scalar_quantity("spon_curvature", dg.readData(
#     outputMesh, 'vertex', 'spon_curvature'), enabled=True)
# ps.show()
