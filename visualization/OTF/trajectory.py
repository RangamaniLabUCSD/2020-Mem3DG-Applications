import netCDF4 as nc
import sys
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pymem3dg as dg
import imp


def matplotlibStyle():
    """ Formatting style of matplotlib """
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


def sizeOf(trajnc):
    with nc.Dataset(trajnc) as ds:
        return ds.groups["Trajectory"].dimensions['frame'].size


def constructSystemByMatrix(parameters, trajnc, frame):
    with nc.Dataset(trajnc) as ds:
        time = np.array(ds.groups['Trajectory'].variables['time'][frame])
        coordinates = np.array(
            ds.groups['Trajectory'].variables['coordinates'][frame])
        topology = np.array(
            ds.groups['Trajectory'].variables['topology'][frame])
        coordinates = np.reshape(coordinates, (-1, 3))
        topology = np.reshape(topology, (-1, 3))
        return dg.System(topology, coordinates, parameters)


def constructSystem(parameters, trajnc, frame):
    return dg.System(trajnc, frame, parameters, 0, True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data", help="traj.nc and parameter setup file used for visualization of data, such as energy trajectory", type=str,  nargs=2)
    args = parser.parse_args()

    """ parse the netCDF trajectory file """
    trajFile = args.data[1]

    """ parse the parameter.py file """
    parameterFile = imp.load_source("module.name", args.data[0])

    """ initialize numpy array """
    frameNum = sizeOf(trajFile)
    time = np.zeros([frameNum, 1])
    kineticEnergy = np.zeros([frameNum, 1])
    potentialEnergy = np.zeros([frameNum, 1])
    externalWork = np.zeros([frameNum, 1])
    totalEnergy = np.zeros([frameNum, 1])

    dirichletEnergy = np.zeros([frameNum, 1])
    bendingEnergy = np.zeros([frameNum, 1])
    surfaceEnergy = np.zeros([frameNum, 1])
    pressureEnergy = np.zeros([frameNum, 1])
    adsorptionEnergy = np.zeros([frameNum, 1])

    """ loop over trajectory and recover data """
    for frame in range(frameNum):
        system = constructSystem(parameterFile.parameters(), trajFile, frame)
        system.computeTotalEnergy()
        time[frame] = system.time
        kineticEnergy[frame] = system.energy.kineticEnergy
        if frame != 0:
            externalWork[frame] = externalWork[frame-1] + \
                system.computeIntegratedPower(time[frame] - time[frame-1])
        potentialEnergy[frame] = system.energy.potentialEnergy
        dirichletEnergy[frame] = system.energy.dirichletEnergy
        bendingEnergy[frame] = system.energy.bendingEnergy
        surfaceEnergy[frame] = system.energy.surfaceEnergy
        pressureEnergy[frame] = system.energy.pressureEnergy
        adsorptionEnergy[frame] = system.energy.adsorptionEnergy
    totalEnergy = potentialEnergy + kineticEnergy - externalWork

    """ plotting """
    fig, axs = plt.subplots(1,2)
    fig.set_size_inches(7, 3)
    # plt.subplots_adjust(left=0.164, bottom=0.07, right=0.988, top=0.988)
    matplotlibStyle()
    axs[0].plot(time, kineticEnergy, label='$E_{kinetic}$')
    axs[0].plot(time, potentialEnergy, label='$E_{potential}$')
    axs[0].plot(time, externalWork, label='$W$')
    axs[0].plot(time, totalEnergy, label='$E_{total}$')
    axs[0].legend()

    axs[1].plot(time, bendingEnergy, label='$E_b$')
    axs[1].plot(time, dirichletEnergy, label='$E_d$')
    axs[1].plot(time, surfaceEnergy, label='$E_s$')
    axs[1].plot(time, pressureEnergy, label='$E_p$')
    axs[1].plot(time, adsorptionEnergy, label='$E_a$')
    axs[1].legend()
    
    plt.tight_layout()
    # plt.savefig("energy.pdf", transparent=True)
    # plt.savefig("energy.png", transparent=True, dpi=1200)
    plt.show()
