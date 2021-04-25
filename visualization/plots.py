import netCDF4 as nc
import sys
import json
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np


def plot(trajnc, figureFile, show=False, save=False):
    # Set unit of data
    kBT = 1.38e-23 * 310 * 1e15  # 1e-15 J
    Pa = 1e-3  # kPa
    pN = 1e-3  # nN

    # Read data from Trajectory file
    ds = nc.Dataset(trajnc)
    issmoothmask = np.array(ds.variables['issmooth'])
    issmoothmask = [i for i, x in enumerate(issmoothmask) if not x]
    # issmoothmask = []
    totalenergy = np.delete(
        np.array(ds.variables['totalenergy'])/kBT, issmoothmask)
    bendenergy =  np.delete(np.array(ds.variables['bendenergy'])/kBT, issmoothmask)
    surfenergy =  np.delete(np.array(ds.variables['surfenergy'])/kBT, issmoothmask)
    pressenergy =  np.delete(np.array(ds.variables['pressenergy'])/kBT, issmoothmask)
    lineenergy =  np.delete(np.array(ds.variables['lineenergy'])/kBT, issmoothmask)
    kineenergy =  np.delete(np.array(ds.variables['kineenergy'])/kBT, issmoothmask)
    chemenergy =  np.delete(np.array(ds.variables['chemenergy'])/kBT, issmoothmask)

    # norm data
    l1errornorm =  np.delete(np.array(ds.variables['l1errornorm'])/Pa, issmoothmask)
    l1bendnorm =  np.delete(np.array(ds.variables['l1bendnorm'])/Pa, issmoothmask)
    l1surfnorm =  np.delete(np.array(ds.variables['l1surfnorm'])/Pa, issmoothmask)
    l1pressnorm =  np.delete(np.array(ds.variables['l1pressnorm'])/Pa, issmoothmask)
    l1linenorm =  np.delete(np.array(ds.variables['l1linenorm'])/Pa, issmoothmask)

    # geometric data
    # surfarea =  np.delete(np.array(ds.variables['surfacearea']), issmoothmask)
    # volume =  np.delete(np.array(ds.variables['volume']), issmoothmask)
    # refsurfarea =  np.delete(np.array(ds.variables['refsurfarea']), issmoothmask)
    # refvolume =  np.delete(np.array(ds.variables['refvolume']), issmoothmask)
    height =  np.delete(np.array(ds.variables['height']), issmoothmask)

    # time
    time =  np.delete(np.array(ds.variables['time']), issmoothmask)

    # Processed data
    # dsurfarea = surfarea / refsurfarea - 1
    # dvolume = volume / refvolume - 1

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
    axs[3].set_xlabel('Time')
    axs[0].set_ylabel('Energy ($k_B T$)')
    axs[1].set_ylabel('Energy ($k_B T$)')
    axs[2].set_ylabel('$L_1$ Norm ($Pa$)')
    axs[3].set_ylabel('Geometry')

    # line graph
    te = axs[0].plot(time, totalenergy, label='$E$')
    ke = axs[0].plot(time, kineenergy, label='$E_{kinetic}$')
    pe = axs[0].plot(time, totalenergy - kineenergy, label='$E_{potential}$')
    axs[0].legend()
    axs[0].set_xticklabels([])

    be = axs[1].plot(time, bendenergy, label='$E_{b}$')
    se = axs[1].plot(time, surfenergy, label='$E_{s}$')
    pse = axs[1].plot(time, pressenergy, label='$E_{p}$')
    ce = axs[1].plot(time, chemenergy, label='$E_{c}$')
    le = axs[1].plot(time, lineenergy, label='$E_{l}$')
    axs[1].legend()
    axs[1].set_xticklabels([])

    l1 = axs[2].plot(time, l1errornorm, label="$e$")
    l1_bend = axs[2].plot(time, l1bendnorm, label="$e_{b}$")
    l1_surf = axs[2].plot(time, l1surfnorm, label="$e_{s}$")
    l1_press = axs[2].plot(time, l1pressnorm, label="$e_{p}$")
    l1_line = axs[2].plot(time, l1linenorm, label="$e_{l}$")
    axs[2].legend()
    axs[2].set_xticklabels([])

    # A = axs[3].plot(time, dsurfarea, label='$A$')
    # V = axs[3].plot(time, dvolume, label='$V$')
    h = axs[3].plot(time, height, label='$h$')
    axs[3].legend()
    # axs[3].ticklabel_format(axis='y', style='sci',
    #                         scilimits=[-3, -3], useMathText=True)

    # plt.tight_layout()

    if save:
        plt.savefig(figureFile, transparent=True)
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
