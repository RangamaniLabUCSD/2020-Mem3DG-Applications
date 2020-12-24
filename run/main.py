import sys
import json
import os
import optparse
import shutil
import pymem3dg
from pathlib import Path

# parse the config file and options
optparser = optparse.OptionParser()
optparser.add_option('-p', '--ply', dest="ply", type="string",
                     help="input mesh and reference mesh files in .ply format")
optparser.add_option('-n', '--nc', dest="nc", type="string",
                     help="input trajectory file in .nc format")
# optparser.add_option('-s', '--swp', dest="swp", type="string",
#                      help="feedforward sweep")
(options, args) = optparser.parse_args()
# argparser = argparse.ArgumentParser()
# parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
# args = argparser.parse_args()
if (options.ply != None):
    configFile = options.ply
elif(options.nc != None):
    configFile = options.nc
# elif(options.swp != None):
#     configFile = options.swp

# read the config file
with open(configFile) as f:
    data = json.load(f)
dep = data["parameters"]["dependencies"]
io = data["parameters"]["I/O"]
opt = data["parameters"]["options"]
var = data["parameters"]["variables"]
prop = data["parameters"]["properties"]
inte = data["parameters"]["integration"]

# configure path for different platforms
io['refMesh'] = os.fspath(os.path.abspath(Path(io['refMesh'])))
io['inputMesh'] = os.fspath(os.path.abspath(Path(io['inputMesh'])))
io['outputDir'] = os.fspath(os.path.abspath(Path(io['outputDir'])))
io['trajFile'] = os.fspath(os.path.abspath(Path(io['trajFile'])))

# make I/O directories if not exist
cwd = os.getcwd()
IDir = os.path.join(cwd, os.path.split(io["inputMesh"])[0])
ODir = os.path.join(cwd, io["outputDir"])
if not os.path.exists(IDir):
    os.mkdir(IDir)
if not os.path.exists(ODir):
    os.mkdir(ODir)

# copy the config.json & viewer to the outputDir
shutil.copyfile(configFile, os.path.join(io["outputDir"], "config.json"))
shutil.copyfile(dep["viewer.py"], os.path.join(io["outputDir"], "viewer.py"))
shutil.copyfile(dep["plots.py"], os.path.join(io["outputDir"], "plots.py"))

# create starting mesh
if io["generateGeometry"] == True:
    pymem3dg.genIcosphere(nSub=io["nSub"], path=io["refMesh"], R=io["R"])
# elif io["inputMesh"] == "UVsphere.ply":
#     pymem3dg.genUVsphere( nSub = io["nSub"]
#                     , path = io["inputMesh"])

# run simulation
if (options.ply != None):
    pymem3dg.driver_ply(verbosity=io["verbosity"],
                        inputMesh=io["inputMesh"],
                        outputDir=io["outputDir"],
                        refMesh=io["refMesh"],
                        nSub=io["nSub"],

                        isTuftedLaplacian=opt["isTuftedLaplacian"],
                        mollifyFactor=opt["mollifyFactor"],
                        isVertexShift=opt["isVertexShift"],
                        isProtein=opt["isProtein"],

                        H0=var["H0*R"] / io["R"],
                        sharpness=var["sharpness"],
                        r_H0=var["r_H0"],
                        Vt=var["Vt"],
                        pt=var["pt"],
                        Kf=var["Kf"],
                        conc=var["conc"],
                        height=var["height"],

                        Kb=prop["Kb"],
                        eta = prop["eta"],
                        Ksg=prop["Ksg"],
                        Kv=prop["Kv"],
                        epsilon=prop["epsilon"],
                        Bc=prop["Bc"],
                        Kse=prop["Kse"],
                        Ksl=prop["Ksl"],
                        Kst=prop["Kst"],
                        kt=prop["kt"],
                        gamma=prop["gamma"],

                        radius=inte["radiusOfIntegration"],
                        h=inte["h"],
                        T=inte["T"],
                        eps=inte["eps"],
                        tSave=inte["tSave"],
                        integration = inte["method"],
                        isBacktrack = inte["options"]["isBacktrack"],
                        rho = inte["options"]["rho"],
                        c1 = inte["options"]["c1"],
                        isAugmentedLagrangian = inte["options"]["isAugmentedLagrangian"])
elif (options.nc != None):
    pymem3dg.driver_nc( verbosity=io["verbosity"],
                        trajFile=io["trajFile"],
                        startingFrame=io["startingFrame"],
                        outputDir=io["outputDir"],

                        isTuftedLaplacian=opt["isTuftedLaplacian"],
                        mollifyFactor=opt["mollifyFactor"],
                        isVertexShift=opt["isVertexShift"],
                        isProtein=opt["isProtein"],

                        H0=var["H0*R"] / io["R"],
                        sharpness=var["sharpness"],
                        r_H0=var["r_H0"],
                        Vt=var["Vt"],
                        pt=var["pt"],
                        Kf=var["Kf"],
                        conc=var["conc"],
                        height=var["height"],

                        Kb=prop["Kb"],
                        eta = prop["eta"],
                        Ksg=prop["Ksg"],
                        Kv=prop["Kv"],
                        epsilon=prop["epsilon"],
                        Bc=prop["Bc"],
                        Kse=prop["Kse"],
                        Ksl=prop["Ksl"],
                        Kst=prop["Kst"],
                        kt=prop["kt"],
                        gamma=prop["gamma"],

                        radius=inte["radiusOfIntegration"],
                        h=inte["h"],
                        T=inte["T"],
                        eps=inte["eps"],
                        tSave=inte["tSave"],
                        integration = inte["method"],
                        isBacktrack = inte["options"]["isBacktrack"],
                        rho = inte["options"]["rho"],
                        c1 = inte["options"]["c1"],
                        isAugmentedLagrangian = inte["options"]["isAugmentedLagrangian"])