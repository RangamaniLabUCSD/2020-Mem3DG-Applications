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
                     help="input mesh and remesh file in .ply format")
optparser.add_option('-n', '--nc', dest="nc", type="string",
                     help="input trajectory file in .nc format")
(options, args) = optparser.parse_args()
# argparser = argparse.ArgumentParser()
# parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
# args = argparser.parse_args()
if (options.ply != None):
    configFile = options.ply
elif(options.nc != None):
    configFile = options.nc

# read the config file
with open(configFile) as f:
    data = json.load(f)
geo = data["parameters"]["geometry"]
prop = data["parameters"]["properties"]
var = data["parameters"]["variables"]
inte = data["parameters"]["integration"]
dep = data["parameters"]["dependencies"]
opt = data["parameters"]["options"]

# configure path for different platforms
geo['refMesh'] = os.fspath(os.path.abspath(Path(geo['refMesh'])))
geo['inputMesh'] = os.fspath(os.path.abspath(Path(geo['inputMesh'])))
geo['outputDir'] = os.fspath(os.path.abspath(Path(geo['outputDir'])))
geo['trajFile'] = os.fspath(os.path.abspath(Path(geo['trajFile'])))

# make I/O directories if not exist
cwd = os.getcwd()
IDir = os.path.join(cwd, os.path.split(geo["inputMesh"])[0])
ODir = os.path.join(cwd, geo["outputDir"])
if not os.path.exists(IDir):
    os.mkdir(IDir)
if not os.path.exists(ODir):
    os.mkdir(ODir)

# copy the config.json & viewer to the outputDir
shutil.copyfile(configFile, os.path.join(geo["outputDir"], "config.json"))
shutil.copyfile(dep["viewer.py"], os.path.join(geo["outputDir"], "viewer.py"))
shutil.copyfile(dep["plots.py"], os.path.join(geo["outputDir"], "plots.py"))

# create starting mesh
if geo["generateGeometry"] == True:
    pymem3dg.genIcosphere(nSub=geo["nSub"], path=geo["refMesh"], R=geo["R"])
# elif geo["inputMesh"] == "UVsphere.ply":
#     pymem3dg.genUVsphere( nSub = geo["nSub"]
#                     , path = geo["inputMesh"])

# run simulation
if (options.ply != None):
    pymem3dg.driver_ply(verbosity=geo["verbosity"],
                        inputMesh=geo["inputMesh"],
                        outputDir=geo["outputDir"],
                        refMesh=geo["refMesh"],
                        nSub=geo["nSub"],
                        radius=geo["radiusOfIntegration"],

                        isTuftedLaplacian=opt["isTuftedLaplacian"],
                        mollifyFactor=opt["mollifyFactor"],
                        isVertexShift=opt["isVertexShift"],
                        isProtein=opt["isProtein"],

                        epsilon=var["epsilon"],
                        Bc=var["Bc"],
                        H0=var["H0*R"] / geo["R"],
                        sharpness=var["sharpness"],
                        r_H0=var["r_H0"],
                        Vt=var["Vt"],
                        pt=var["pt"],
                        Kf=var["Kf"],
                        conc=var["conc"],
                        height=var["height"],

                        Kb=prop["Kb"],
                        eta = prop["eta"],
                        Kse=prop["Kse"],
                        Ksl=prop["Ksl"],
                        Kst=prop["Kst"],
                        Ksg=prop["Ksg"],
                        Kv=prop["Kv"],
                        kt=prop["kt"],
                        gamma=prop["gamma"],

                        h=inte["h"],
                        T=inte["T"],
                        eps=inte["eps"],
                        closeZone=inte["closeZone"],
                        increment=inte["increment"],
                        tSave=inte["tSave"],
                        tMollify=inte["tMollify"],
                        errorJumpLim=inte["errorJumpLim"],
                        integration = inte["method"])

elif (options.nc != None):
    pymem3dg.driver_nc(verbosity=geo["verbosity"],
                       trajFile=geo["trajFile"],
                       startingFrame=geo["startingFrame"],
                       outputDir=geo["outputDir"],
                       radius=geo["radiusOfIntegration"],

                       isTuftedLaplacian=opt["isTuftedLaplacian"],
                       mollifyFactor=opt["mollifyFactor"],
                       isVertexShift=opt["isVertexShift"],
                       isProtein=opt["isProtein"],

                       epsilon=var["epsilon"],
                       Bc=var["Bc"],
                       H0=var["H0*R"] / geo["R"],
                       sharpness=var["sharpness"],
                       r_H0=var["r_H0"],
                       Vt=var["Vt"],
                       pt=var["pt"],
                       Kf=var["Kf"],
                       conc=var["conc"],
                       height=var["height"],

                       Kb=prop["Kb"],
                       eta = prop["eta"],
                       Kse=prop["Kse"],
                       Ksl=prop["Ksl"],
                       Kst=prop["Kst"],
                       Ksg=prop["Ksg"],
                       Kv=prop["Kv"],
                       kt=prop["kt"],
                       gamma=prop["gamma"],

                       h=inte["h"],
                       T=inte["T"],
                       eps=inte["eps"],
                       closeZone=inte["closeZone"],
                       increment=inte["increment"],
                       tSave=inte["tSave"],
                       tMollify=inte["tMollify"],
                       errorJumpLim=inte["errorJumpLim"],
                       integration = inte["method"])
