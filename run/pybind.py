import pymem3dg

Kb = 8.22e-5
H0 = 0
sharpness = 0
r_H0 = [-1, -1]
Ksg = 1
Kst = 7
Ksl = 0
Kse = 0
epsilon = -1
Bc = -1
gamma = 0
Kv = 0.1
Vt = -1
Pam = 0.3
Kf = 0
conc = -1
height = 0
radius = 100000
temp = 0
eta = 0
pt = [0, 0, 0]

h = 0.0001
T = 10
eps = 1e-6
tSave = 1
isBacktrack = True
isAdaptiveStep = True

sigma = (2 * gamma * 1.380649e-8 * temp / h)**0.5

isProtein = False
isVertexShift = False
isReducedVolume = False
isLocalCurvature = False
rho = 0.5
c1 = 0
p = pymem3dg.Parameters(Kb,    H0,  sharpness, r_H0, Ksg,    Kst,   Ksl, Kse,
                        Kv,    eta, epsilon,   Bc,   gamma,  Vt,    Pam, temp,
                        sigma, pt,  Kf,        conc, height, radius, 0, 0)
inputMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/slightlyOblate.ply"
refMesh = "/home/cuzhu/2020-Mem3DG-Applications/run/input-file/slightlyOblate.ply"
outputDir = "/home/cuzhu/2020-Mem3DG-Applications/results/vesicle/testAPI"
f = pymem3dg.System(inputMesh, refMesh, 1, p, isReducedVolume,
                    isProtein, isLocalCurvature, isVertexShift)
f.computeFreeEnergy()
# print(f.E.totalE)
# print(f.computeChemicalPotential())
# print(f.computeBendingPressure())
# print(f.P.r_H0)
f.computeBendingPressure()
print(f.getBendingPressure())
print(f.getProteinDensity())
print(f.getVertexPositionMatrix())
pymem3dg.euler(f, h, T, tSave, 0, 3, outputDir, isBacktrack, rho, c1, isAdaptiveStep)
# pymem3dg.velocityVerlet(f, h, T, tSave, 0, 3, isAdaptiveStep,outputDir)
