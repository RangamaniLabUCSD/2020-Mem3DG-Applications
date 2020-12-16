import netCDF4 as nc
import sys
import json
import argparse
import os
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import pymem3dg
import threading

# parse the config file
parser = argparse.ArgumentParser()
# parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
parser.add_argument(
    "data", help="traj.nc file used for visualization of data, such as energy trajectory", type=str)
args = parser.parse_args()
ext = os.path.splitext(args.data)[1][1:]

ds = nc.Dataset(args.data)
bendenergy = ds.variables['bendenergy']
surfenergy = ds.variables['surfenergy']
pressenergy = ds.variables['pressenergy']
kineenergy = ds.variables['kineenergy']
chemenergy = ds.variables['chemenergy']
lineenergy = ds.variables['lineenergy']
totalenergy = ds.variables['totalenergy']

l2errornorm = ds.variables['l2errornorm']

time = ds.variables['time']

fig, axs = plt.subplots(3)
fig.suptitle('Energy trajectory')
axs[2].set_xlabel('time (s)')
axs[0].set_ylabel('energy (10^(-15) J)')
axs[1].set_ylabel('energy (10^(-15) J)')
axs[2].set_ylabel('L2 Error Norm (10^(-9) N)')

te = axs[0].plot(time, totalenergy, label='total')
ke = axs[0].plot(time, kineenergy, label='kinetic')
pe = axs[0].plot(time, np.array(totalenergy) -
                  np.array(kineenergy), label='potential')
axs[0].legend()

be = axs[1].plot(time, bendenergy, label='bending')
se = axs[1].plot(time, surfenergy, label='surface')
pse = axs[1].plot(time, pressenergy, label='pressure')
ce = axs[1].plot(time, chemenergy, label='chemical')
le = axs[1].plot(time, lineenergy, label='line')
axs[1].legend()

l2 = axs[2].plot(time, l2errornorm, label = 'L2 error')
axs[2].legend()

plt.show()

# fig = plt.figure()
# plt.plot( 'x', 'y1', data=df, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
# plt.plot( 'x', 'y2', data=df, marker='', color='olive', linewidth=2)
# plt.plot( 'x', 'y3', data=df, marker='', color='olive', linewidth=2, linestyle='dashed', label="toto")

#plt.plot(dt_time[time_idx], air[time_idx, lat_idx, lon_idx], c='b', marker='o')
# ax.legend((te, ke, pe), ('total', 'kinetic', 'potential'), loc='upper right', shadow=False)
# axs[1].legend([be, se, pse, ce], ['bending', 'surface', 'pressure', 'chemical'], loc='upper right', shadow=False)
