import sys
import json
import os
import optparse
import shutil
import pymem3dg
import numpy as np
from pathlib import Path
from main import configParse, genPlots

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
                                     isAugmentedLagrangian=inte["options"]["isAugmentedLagrangian"])


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
                                    isAugmentedLagrangian=inte["options"]["isAugmentedLagrangian"])

if __name__ == "__main__":
    # parse the command line option
    optparser = optparse.OptionParser()
    optparser.add_option('-p', '--ply', dest="ply", type="string",
                         help="input mesh and reference mesh files in .ply format")
    optparser.add_option('-n', '--nc', dest="nc", type="string",
                         help="input trajectory file in .nc format")
    (options, args) = optparser.parse_args()
    # argparser = argparse.ArgumentParser()
    # parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
    # args = argparser.parse_args()

    # parse the config file
    dep, io, opt, var, prop, inte = configParse(options, args)

    # run simulation
    if (options.ply != None):
        plySweep(dep, io, opt, var, prop, inte)
    elif(options.nc != None):
        ncSweep(dep, io, opt, var, prop, inte)

    # generate plots based on netcdf trajectory
    genPlots(io)
