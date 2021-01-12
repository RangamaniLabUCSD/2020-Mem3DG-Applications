import sys
import json
import os
import argparse
import shutil
import pymem3dg
import numpy as np
from pathlib import Path
from main import configParse, genPlots

import multiprocessing
import subprocess
from tqdm.contrib.concurrent import process_map

import concurrent.futures

def plySweep(dep, io, opt, var, prop, inte):
    '''function run the forward sweep function starting with .ply mesh'''
    # create starting mesh
    if io["generateGeometry"] == True:
        pymem3dg.genIcosphere(nSub=io["nSub"], path=io["refMesh"], R=io["R"])
    # elif io["inputMesh"] == "UVsphere.ply":
    #     pymem3dg.genUVsphere( nSub = io["nSub"]
    #                     , path = io["inputMesh"])
    return pymem3dg.forwardsweep_ply(inputMesh=io["inputMesh"],
                                     outputDir=io["outputDir"],
                                     refMesh=io["refMesh"],
                                     nSub=io["nSub"],

                                     isVertexShift=opt["isVertexShift"],
                                     isProtein=opt["isProtein"],
                                     isLocalCurvature=opt["isLocalCurvature"],
                                     isReducedVolume=opt["isReducedVolume"],

                                     H0=np.array(var["H0*R"]) / io["R"],
                                     sharpness=var["sharpness"],
                                     r_H0=var["r_H0"],
                                     Vt=var["Vt"],
                                     cam=var["cam"],
                                     pt=var["pt"],
                                     Kf=var["Kf"],
                                     conc=var["conc"],
                                     height=var["height"],

                                     Kb=prop["Kb"],
                                     eta=prop["eta"],
                                     Ksg=prop["Ksg"],
                                     Kv=prop["Kv"],
                                     epsilon=prop["epsilon"],
                                     Bc=prop["Bc"],
                                     Kse=prop["Kse"],
                                     Ksl=prop["Ksl"],
                                     Kst=prop["Kst"],
                                     temp=prop["temp"],
                                     gamma=prop["gamma"],

                                     radius=inte["radiusOfIntegration"],
                                     h=inte["h"],
                                     T=inte["T"],
                                     eps=inte["eps"],
                                     tSave=inte["tSave"],
                                     isBacktrack=inte["options"]["isBacktrack"],
                                     rho=inte["options"]["rho"],
                                     c1=inte["options"]["c1"],
                                     ctol=inte["options"]["ctol"],
                                     isAugmentedLagrangian=inte["options"]["isAugmentedLagrangian"],
                                     isAdaptiveStep=inte["options"]["isAdaptiveStep"])


def ncSweep(dep, io, opt, var, prop, inte):
    '''function that runs the forward sweep function starting with data stored in netcdf file'''
    return pymem3dg.forwardsweep_nc(trajFile=io["trajFile"],
                                    startingFrame=io["startingFrame"],
                                    outputDir=io["outputDir"],

                                    isVertexShift=opt["isVertexShift"],
                                    isProtein=opt["isProtein"],
                                    isLocalCurvature=opt["isLocalCurvature"],
                                    isReducedVolume=opt["isReducedVolume"],

                                    H0=np.array(var["H0*R"]) / io["R"],
                                    sharpness=var["sharpness"],
                                    r_H0=var["r_H0"],
                                    Vt=var["Vt"],
                                    cam=var["cam"],
                                    pt=var["pt"],
                                    Kf=var["Kf"],
                                    conc=var["conc"],
                                    height=var["height"],

                                    Kb=prop["Kb"],
                                    eta=prop["eta"],
                                    Ksg=prop["Ksg"],
                                    Kv=prop["Kv"],
                                    epsilon=prop["epsilon"],
                                    Bc=prop["Bc"],
                                    Kse=prop["Kse"],
                                    Ksl=prop["Ksl"],
                                    Kst=prop["Kst"],
                                    temp=prop["temp"],
                                    gamma=prop["gamma"],

                                    radius=inte["radiusOfIntegration"],
                                    h=inte["h"],
                                    T=inte["T"],
                                    eps=inte["eps"],
                                    tSave=inte["tSave"],
                                    isBacktrack=inte["options"]["isBacktrack"],
                                    rho=inte["options"]["rho"],
                                    c1=inte["options"]["c1"],
                                    ctol=inte["options"]["ctol"],
                                    isAugmentedLagrangian=inte["options"]["isAugmentedLagrangian"],
                                    isAdaptiveStep=inte["options"]["isAdaptiveStep"])


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
