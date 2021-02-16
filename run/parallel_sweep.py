import sys
import json
import os
import optparse
import shutil
import numpy as np
from pathlib import Path

import pymem3dg
from main import configParse, genPlots, plyRun, ncRun

import multiprocessing
import subprocess
from tqdm.contrib.concurrent import process_map

import concurrent.futures
def worker(args):
  # read all simulation parameters
    (dep, io, opt, var, prop, inte) = args[0]

    # read the thread specific sweeping values and assign to parameters
    H0 = args[1]
    cam = args[2]
    Vt = args[3]
    Ksg = args[4]
    var["H0"] = H0
    var["cam"] = cam
    var["Vt"] = Vt
    prop["Ksg"] = Ksg

    # construct the subfolder for each run
    subFolder = args[5]
    subFolderPath = os.path.join(io["outputDir"], subFolder)
    io["outputDir"] = subFolderPath
    if not os.path.exists(subFolderPath):
        os.mkdir(subFolderPath)

    # run the sweep
    # plyRun(dep, io, opt, var, prop, inte)

    # run simulation
    if (options.ply != None):
        plyRun(dep, io, opt, var, prop, inte)
    elif(options.nc != None):
        ncRun(dep, io, opt, var, prop, inte)


if __name__ == "__main__":
    
    # # parse the json argument
    # parser = argparse.ArgumentParser()
    # parser.add_argument(
    #     "config", help="configuration file (.json) used for the simulation", type=str)
    # configFile = parser.parse_args().config

    # parse the command line option
    optparser = optparse.OptionParser()
    optparser.add_option('-p', '--ply', dest="ply", type="string",
                         help="input mesh and reference mesh files in .ply format")
    optparser.add_option('-n', '--nc', dest="nc", type="string",
                         help="input trajectory file in .nc format")
    (options, args) = optparser.parse_args()
    if (options.ply != None):
        configFile = options.ply
    elif(options.nc != None):
        configFile = options.nc

    # parse the config file
    dep, io, opt, var, prop, inte = configParse(configFile)

    # parallel run the axis of cam for each H0[i]
    jobs = []
    for H0 in np.array(var["H0"]):
        for cam in np.array(var["cam"]):
            for Vt in np.array(var["Vt"]):
                for Ksg in np.array(prop["Ksg"]):
                    subFolder = f'H_{H0}_c_{cam}_V_{Vt}_Ksg_{Ksg}/'
                    jobs.append(((dep, io, opt, var, prop, inte),
                                H0, cam, Vt, Ksg, subFolder))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(worker, jobs)

    # generate plots based on netcdf trajectory
    # genPlots(io)
