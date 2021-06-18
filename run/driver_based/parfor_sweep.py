import sys
import json
import os
import argparse
import shutil
import numpy as np
from pathlib import Path

import pymem3dg
from main import configParse, genPlots
from forward_sweep import plySweep, ncSweep

import multiprocessing
import subprocess
from tqdm.contrib.concurrent import process_map

import concurrent.futures

def worker(args):
    # read all simulation parameters
    (dep, io, opt, var, prop, inte) = args[0]

    # read the thread specific H0 value and assign to parameters
    H0 = args[1]
    var["H0*R"] = [H0]

    # read starting Cam used to identify starting trajFile 
    cam_ = args[2]
    startingCam = cam_[0]
    trajFileName = f'traj_H_{int(H0 * 100)}_VP_{int(startingCam * 100)}.nc'
    io["trajFile"] = os.path.join(io["outputDir"],trajFileName)

    # run the rest in the axis of cam 
    var["cam"] = cam_[1:]

    # the starting frame is the last frame of the trajFile
    io["startingFrame"] = -1

    # run the sweep
    ncSweep(dep, io, opt, var, prop, inte)

if __name__ == "__main__":
    # parse the json argument
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
    configFile = parser.parse_args().config

    # parse the config file
    dep, io, opt, var, prop, inte = configParse(configFile)

    # save vector of sweeping cam 
    cam_ = np.array(var["cam"])

    # first run the axis of H0 with cam[0]
    var["cam"] = [cam_[0]]
    plySweep(dep, io, opt, var, prop, inte)

    # parallel run the axis of cam for each H0[i]
    var["cam"] = cam_
    jobs = []
    for H0 in np.array(var["H0*R"]):        
        jobs.append(((dep, io, opt, var, prop, inte), H0, cam_))
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(worker, jobs)
    
    # generate plots based on netcdf trajectory
    # genPlots(io)
