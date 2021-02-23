import pymem3dg
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/slightlyOblate.ply"
# refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/slightlyOblate.ply"
# outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/vesicle/testAPI"
inputMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//slightlyOblate.ply"
refMesh = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//run//input-file//slightlyOblate.ply"
outputDir = "C://Users//Kieran//Dev//2020-Mem3DG-Applications//results//vesicle//testConvergence"

isProtein = False
isVertexShift = False
isReducedVolume = False
isLocalCurvature = True

nSub = 0

p = pymem3dg.Parameters()
p.Kb = 8.22e-5
p.H0 = 3
p.sharpness = 10000000
p.r_H0 = [1, 1]
p.eta = 0.0005
p.Ksg = 0.05
p.Kst = 0.1
p.Ksl = 0
p.Kse = 0
p.epsilon = -1
p.Bc = -1
p.Kv = 0.1
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


def myExpFunc(x, a, b):
    return a * np.power(x, b)


# Set unit of data
kBT = 1.38e-23 * 310 * 1e15  # 1e-15 J
pJ = 1e3  # fJ = 1e-15 J
Pa = 1e-3  # kPa

nSubSize = 5
nvertices = np.zeros(nSubSize)
l2bendnorm = np.zeros(nSubSize)
l2surfnorm = np.zeros(nSubSize)
l2pressnorm = np.zeros(nSubSize)
l2linenorm = np.zeros(nSubSize)

bendenergy = np.zeros(nSubSize)
surfenergy = np.zeros(nSubSize)
pressenergy = np.zeros(nSubSize)
lineenergy = np.zeros(nSubSize)

for i in range(nSubSize):
    n = i+1
    faceMatrix, vertexMatrix = pymem3dg.getIcosphere(n, 1.3)
    faceMatrix, refVertexMatrix = pymem3dg.getIcosphere(n, 1)
    g = pymem3dg.System(faceMatrix, vertexMatrix, refVertexMatrix, nSub, p, isReducedVolume,
                        isProtein, isLocalCurvature, isVertexShift)
    g.computeFreeEnergy()
    g.computeAllForces()
    nvertices[i] = len(np.array(g.getVertexPositionMatrix()))
    l2bendnorm[i] = g.computeL2Norm(
        g.getLumpedMassMatrix() * g.getBendingPressure())
    l2surfnorm[i] = g.computeL2Norm(
        g.getLumpedMassMatrix() * g.getCapillaryPressure())
    l2pressnorm[i] = g.computeL2Norm(
        g.getLumpedMassMatrix() * g.getInsidePressure())
    l2linenorm[i] = g.computeL2Norm(
        g.getLumpedMassMatrix() * g.getLineTensionPressure())

    bendenergy[i] = g.E.BE
    surfenergy[i] = g.E.sE
    pressenergy[i] = g.E.pE
    lineenergy[i] = g.E.lE


print(g.getLineTension())
print(bendenergy)
print(surfenergy)
print(pressenergy)
print(lineenergy)
print(abs(l2bendnorm - l2bendnorm[nSubSize - 1])
      [:nSubSize - 1])
# l2bendnorm = l2bendnorm / Pa
# l2surfnorm = l2surfnorm / Pa
# l2pressnorm = l2pressnorm / Pa

# Visualization preference
# Font:
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
#mpl.rcParams.update({'font.size': 8})
SMALL_SIZE = 13
MEDIUM_SIZE = 18
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('pdf', fonttype=42)
# plt.subplots_adjust(left=0.164, bottom=0.07, right=0.988, top=0.988)


ax = plt.subplot(1, 1, 1)
# ax.plot(1/nvertices[:nSubSize-1], abs(l2bendnorm - l2bendnorm[nSubSize-1])[:nSubSize-1], label="bending")
# ax.plot(1/nvertices[:nSubSize-1], abs(l2surfnorm - l2surfnorm[nSubSize-1])[:nSubSize-1], label="capillary")
# ax.plot(1/nvertices[:nSubSize-1], abs(l2pressnorm - l2pressnorm[nSubSize-1])[:nSubSize-1], label="pressure")
newX = np.logspace(-1.6, -1.1, base=10)

popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(l2bendnorm - l2bendnorm[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(l2bendnorm - l2bendnorm[nSubSize-1])[:nSubSize-1],
          '-o',  label="$f_b, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(l2surfnorm - l2surfnorm[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(l2surfnorm - l2surfnorm[nSubSize-1])[:nSubSize-1],
          '-o', label="$f_s, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(l2pressnorm - l2pressnorm[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(l2pressnorm - l2pressnorm[nSubSize-1])[:nSubSize-1],
          '-o', label="$f_p, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(l2linenorm - l2linenorm[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(l2linenorm - l2linenorm[nSubSize-1])[:nSubSize-1],
          '-o', label="$f_l, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

ax.legend()
ax.set_ylabel("$e_f$(kPa)")
ax.set_xlabel("$L$ (# Vertices$^{-1/2}$)")
plt.tight_layout()
plt.savefig("pressure_conv.png")
ax.clear()

newX = np.logspace(-2.5, -1.3, base=10)
popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(bendenergy - bendenergy[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(bendenergy - bendenergy[nSubSize-1])[:nSubSize-1],
          '-o', label="$E_b, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(surfenergy - surfenergy[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(surfenergy - surfenergy[nSubSize-1])[:nSubSize-1],
          '-o', label="$E_s, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(pressenergy - pressenergy[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(pressenergy - pressenergy[nSubSize-1])[:nSubSize-1],
          '-o', label="$E_p, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(lineenergy - lineenergy[nSubSize-1]))
ax.loglog(np.power(nvertices, -0.5)[:nSubSize-1],
          abs(lineenergy - lineenergy[nSubSize-1])[:nSubSize-1],
          '-o', label="$E_l, {0:.1f}L^{{{1:.1f}}}$".format(*popt))


ax.legend()
ax.set_ylabel("$e_E(k_BT)$")
ax.set_xlabel("$L$ (# Vertices$^{-1/2}$)")
plt.tight_layout()
plt.savefig("energy_conv.png")



# plt.show()
