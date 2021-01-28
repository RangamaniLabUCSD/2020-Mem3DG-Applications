import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def myExpFunc(x, a, b):
    return a * np.power(x, b)


# Set unit of data
kBT = 1.38e-23 * 310 * 1e15  # 1e-15 J
pJ = 1e3 # fJ = 1e-15 J
Pa = 1e-3  # kPa

nvertices = np.zeros(6)
l2bendnorm = np.zeros(6)
l2surfnorm = np.zeros(6)
l2pressnorm = np.zeros(6)

bendenergy = np.zeros(6)
surfenergy = np.zeros(6)
pressenergy = np.zeros(6)


for i in range(6):
    ds = nc.Dataset("traj_{}.nc".format(i))
    nvertices[i] = len(np.array(ds.variables['refcoordinates']))
    l2bendnorm[i] = np.array(ds.variables['l2bendnorm'])[0]
    l2surfnorm[i] = np.array(ds.variables['l2surfnorm'])[0]
    l2pressnorm[i] = np.array(ds.variables['l2pressnorm'])[0]
    bendenergy[i] = np.array(ds.variables['bendenergy'])[0]/kBT
    surfenergy[i] = np.array(ds.variables['surfenergy'])[0]/kBT
    pressenergy[i] = np.array(ds.variables['pressenergy'])[0]/kBT

print(l2bendnorm)
print(l2surfnorm)
print(l2pressnorm)
print(abs(l2bendnorm - l2bendnorm[5])
          [:5])
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
# ax.plot(1/nvertices[:5], abs(l2bendnorm - l2bendnorm[5])[:5], label="bending")
# ax.plot(1/nvertices[:5], abs(l2surfnorm - l2surfnorm[5])[:5], label="capillary")
# ax.plot(1/nvertices[:5], abs(l2pressnorm - l2pressnorm[5])[:5], label="pressure")
ax.loglog(np.power(nvertices, -0.5)[:5], abs(l2bendnorm - l2bendnorm[5])
          [:5], '-o', label="$f_b$")
ax.loglog(np.power(nvertices, -0.5)[:5], abs(l2surfnorm - l2surfnorm[5])
          [:5], '-o', label="$f_s$")
ax.loglog(np.power(nvertices, -0.5)[:5], abs(l2pressnorm -
                                             l2pressnorm[5])[:5], '-o', label="$f_p$")

newX = np.logspace(-2.5, -1.3, base=10)
popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(l2bendnorm - l2bendnorm[5]))
plt.plot(newX, myExpFunc(newX, *popt), '--',
         label="$f_b, {0:.1f}L^{{{1:.1f}}}$".format(*popt))
popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(l2surfnorm - l2surfnorm[5]))
plt.plot(newX, myExpFunc(newX, *popt), '--',
         label="$f_s, {0:.1f}L^{{{1:.1f}}}$".format(*popt))
popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(l2pressnorm - l2pressnorm[5]))
plt.plot(newX, myExpFunc(newX, *popt), '--',
         label="$f_p, {0:.1f}L^{{{1:.1f}}}$".format(*popt))

ax.legend()
ax.set_ylabel("$e_f$(kPa)")
ax.set_xlabel("$L$ (# Vertices$^{-1/2}$)")
plt.tight_layout()
plt.savefig("pressure_conv.png")
ax.clear()


ax.loglog(np.power(nvertices, -0.5)[:5], abs(bendenergy - bendenergy[5])
          [:5], '-o', label="$E_b$")
ax.loglog(np.power(nvertices, -0.5)[:5], abs(surfenergy - surfenergy[5])
          [:5], '-o', label="$E_s$")
ax.loglog(np.power(nvertices, -0.5)[:5], abs(pressenergy -
                                             pressenergy[5])[:5], '-o', label="$E_p$")

newX = np.logspace(-2.5, -1.3, base=10)
popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(bendenergy - bendenergy[5]))
plt.plot(newX, myExpFunc(newX, *popt), '--',
         label="$E_b, {0:.1f}L^{{{1:.1f}}}$".format(*popt))
popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(surfenergy - surfenergy[5]))
plt.plot(newX, myExpFunc(newX, *popt), '--',
         label="$E_s, {0:.1f}L^{{{1:.1f}}}$".format(*popt))
popt, pcov = curve_fit(myExpFunc, np.power(
    nvertices, -0.5), abs(pressenergy - pressenergy[5]))
plt.plot(newX, myExpFunc(newX, *popt), '--',
         label="$E_p, {0:.1f}L^{{{1:.1f}}}$".format(*popt))
ax.legend()
ax.set_ylabel("$e_E(k_BT)$")
ax.set_xlabel("$L$ (# Vertices$^{-1/2}$)")
plt.tight_layout()
plt.savefig("energy_conv.png")
plt.show()
