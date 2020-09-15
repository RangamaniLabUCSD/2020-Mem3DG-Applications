import sys
import json
import os
import optparse
import shutil
import pymem3dg

# parse the config file and options
optparser = optparse.OptionParser()
optparser.add_option('-p','--ply', dest = "ply", type = "string", help = "input mesh and remesh file in .ply format")
optparser.add_option('-n', '--nc', dest = "nc", type = "string", help = "input trajectory file in .nc format")
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

# make I/O directories if not exist
cwd = os.getcwd()
IDir = os.path.join(cwd, os.path.split(geo["inputMesh"])[0])
ODir = os.path.join(cwd, geo["outputDir"])
if not os.path.exists(IDir):
    os.mkdir(IDir)
if not os.path.exists(ODir):
    os.mkdir(ODir)

# copy the config.json & viewer to the outputDir
shutil.copyfile(configFile, geo["outputDir"] + "config.json" )
shutil.copyfile("viewer.py", geo["outputDir"] + "viewer.py")

# create starting mesh
if geo["generateGeometry"] == True:
    pymem3dg.genIcosphere( nSub = geo["nSub"] 
                        , path = geo["refMesh"], R = geo["R"])
# elif geo["inputMesh"] == "UVsphere.ply":
#     pymem3dg.genUVsphere( nSub = geo["nSub"] 
#                     , path = geo["inputMesh"])

# run simulation
if (options.ply != None):
    pymem3dg.driver_ply(inputMesh = geo["inputMesh"],
                outputDir = geo["outputDir"],
                refMesh = geo["refMesh"],
                radius = geo["radiusOfIntegration"],

                isTuftedLaplacian = opt["isTuftedLaplacian"],
                mollifyFactor = opt["mollifyFactor"],
                isVertexShift = opt["isVertexShift"],
                isProtein = opt["isProtein"],

                epsilon = var["epsilon"],
                Bc = var["Bc"],
                H0 = var["H0*R"] / geo["R"],
                sharpness = var["sharpness"],
                r_H0 = var["r_H0"],
                Vt = var["Vt"],
                ptInd = var["ptInd"],  
                Kf = var["Kf"],
                conc = var["conc"],
                height = var["height"],
                
                Kb  = prop["Kb"],
                Kse = prop["Kse"],     
                Ksl = prop["Ksl"],
                Kst = prop["Kst"],		
                Ksg = prop["Ksg"],		
                Kv  = prop["Kv"], 	
                kt = prop["kt"], 
                gamma = prop["gamma"], 
                
                h = inte["h"],
                T = inte["T"],
                eps = inte["eps"],
                closeZone = inte["closeZone"],
                increment = inte["increment"],
                tSave = inte["tSave"],
                tMollify = inte["tMollify"])

elif (options.nc != None):
    pymem3dg.driver_nc(trajFile = geo["trajFile"],
                startingFrame = geo["startingFrame"],
                outputDir = geo["outputDir"],
                radius = geo["radiusOfIntegration"],

                isTuftedLaplacian = opt["isTuftedLaplacian"],
                mollifyFactor = opt["mollifyFactor"],
                isVertexShift = opt["isVertexShift"],
                isProtein = opt["isProtein"],

                epsilon = var["epsilon"],
                Bc = var["Bc"],
                H0 = var["H0*R"] / geo["R"],
                sharpness = var["sharpness"],
                r_H0 = var["r_H0"],
                Vt = var["Vt"],
                ptInd = var["ptInd"],  
                Kf = var["Kf"],
                conc = var["conc"],
                height = var["height"],
                
                Kb  = prop["Kb"],
                Kse = prop["Kse"],     
                Ksl = prop["Ksl"],
                Kst = prop["Kst"],		
                Ksg = prop["Ksg"],		
                Kv  = prop["Kv"], 	
                kt = prop["kt"], 
                gamma = prop["gamma"], 
                
                h = inte["h"],
                T = inte["T"],
                eps = inte["eps"],
                closeZone = inte["closeZone"],
                increment = inte["increment"],
                tSave = inte["tSave"],
                tMollify = inte["tMollify"])