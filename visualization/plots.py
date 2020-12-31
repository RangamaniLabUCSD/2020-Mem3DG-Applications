import netCDF4 as nc
import sys
import json
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np


def plot(trajnc, figureFile, show=False, save=False):
    # Read data from Trajectory file
    ds = nc.Dataset(trajnc)
    bendenergy = ds.variables['bendenergy']
    surfenergy = ds.variables['surfenergy']
    pressenergy = ds.variables['pressenergy']
    kineenergy = ds.variables['kineenergy']
    chemenergy = ds.variables['chemenergy']
    lineenergy = ds.variables['lineenergy']
    totalenergy = ds.variables['totalenergy']
    l2errornorm = ds.variables['l2errornorm']
    time = ds.variables['time']

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

    # Figure:
    fig, axs = plt.subplots(3)
    fig.set_size_inches(8, 10)
    plt.subplots_adjust(left=0.164, bottom=0.07, right=0.988, top=0.988)

    # Plotting
    # Title:
    # fig.suptitle('Energy trajectory')

    # Axes label:
    axs[2].set_xlabel('Time ($s$)')
    axs[0].set_ylabel('Energy ($10^{-15} ~J$)')
    axs[1].set_ylabel('Energy ($10^{-15} ~J$)')
    axs[2].set_ylabel('$L_2$ Residual ($10^{-9} ~N$)')

    # line graph
    te = axs[0].plot(time, totalenergy, label='Total')
    ke = axs[0].plot(time, kineenergy, label='Kinetic')
    pe = axs[0].plot(time, np.array(totalenergy) -
                     np.array(kineenergy), label='Potential')
    axs[0].legend()

    be = axs[1].plot(time, bendenergy, label='Bending')
    se = axs[1].plot(time, surfenergy, label='Surface')
    pse = axs[1].plot(time, pressenergy, label='Pressure')
    ce = axs[1].plot(time, chemenergy, label='Chemical')
    le = axs[1].plot(time, lineenergy, label='Line')
    axs[1].legend()

    l2 = axs[2].plot(time, l2errornorm)
    # axs[2].legend()
    if save:
        plt.savefig(figureFile)
    if show:
        plt.show()

    # Archieved
    # fig = plt.figure()
    # plt.plot( 'x', 'y1', data=df, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
    # plt.plot( 'x', 'y2', data=df, marker='', color='olive', linewidth=2)
    # plt.plot( 'x', 'y3', data=df, marker='', color='olive', linewidth=2, linestyle='dashed', label="toto")

    #plt.plot(dt_time[time_idx], air[time_idx, lat_idx, lon_idx], c='b', marker='o')
    # ax.legend((te, ke, pe), ('total', 'kinetic', 'potential'), loc='upper right', shadow=False)
    # axs[1].legend([be, se, pse, ce], ['bending', 'surface', 'pressure', 'chemical'], loc='upper right', shadow=False)


if __name__ == "__main__":
    # Parse the trajectory file
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--save', action='store_true',
                        help="save pdf")
    parser.add_argument(
        "data", help="traj.nc file used for visualization of data, such as energy trajectory", type=str)
    args = parser.parse_args()

    # Run plot()
    plot(trajnc=args.data, figureFile='traj_plot.pdf', show=True, save=args.save)
