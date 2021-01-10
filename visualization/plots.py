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
    bendenergy = np.array(ds.variables['bendenergy'])
    surfenergy = np.array(ds.variables['surfenergy'])
    pressenergy = np.array(ds.variables['pressenergy'])
    kineenergy = np.array(ds.variables['kineenergy'])
    chemenergy = np.array(ds.variables['chemenergy'])
    lineenergy = np.array(ds.variables['lineenergy'])
    totalenergy = np.array(ds.variables['totalenergy'])
    l2errornorm = np.array(ds.variables['l2errornorm'])
    surfarea = np.array(ds.variables['surfacearea'])
    volume = np.array(ds.variables['volume'])
    time = np.array(ds.variables['time'])

    # Processed data 
    surfarea_ = surfarea / surfarea[0]
    refVolume = (surfarea[0] / 4 / np.pi )**(1.5) * 4.0 / 3.0 * np.pi
    volume_ = volume / refVolume

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
    fig, axs = plt.subplots(4)
    fig.set_size_inches(7, 10)
    plt.subplots_adjust(left=0.164, bottom=0.07, right=0.988, top=0.988)

    # Plotting
    # Title:
    # fig.suptitle('Energy trajectory')

    # Axes label:
    axs[2].set_xlabel('Time ($s$)')
    axs[0].set_ylabel('Energy ($10^{-15} ~J$)')
    axs[1].set_ylabel('Energy ($10^{-15} ~J$)')
    axs[2].set_ylabel('$L_2$ Residual ($10^{-9} ~N$)')
    axs[3].set_ylabel('Geometry')

    # line graph
    te = axs[0].plot(time, totalenergy, label='$E$')
    ke = axs[0].plot(time, kineenergy, label='$E_{kinetic}$')
    pe = axs[0].plot(time, totalenergy - kineenergy, label='$E_{potential}$')
    axs[0].legend()

    be = axs[1].plot(time, bendenergy, label='$E_{b}$')
    se = axs[1].plot(time, surfenergy, label='$E_{s}$')
    pse = axs[1].plot(time, pressenergy, label='$E_{p}$')
    ce = axs[1].plot(time, chemenergy, label='$E_{c}$')
    le = axs[1].plot(time, lineenergy, label='$E_{l}$')
    axs[1].legend()

    l2 = axs[2].plot(time, l2errornorm)
    # axs[2].legend()

    A = axs[3].plot(time, surfarea_, label='$A$')
    V = axs[3].plot(time, volume_, label='$V$')
    axs[3].legend()

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
