import sys
import json
import os
import optparse
import shutil
import pymem3dg
from pathlib import Path


def configParse(options, args):
    '''function that parse the .json config file'''
    if (options.ply != None):
        configFile = options.ply
    elif(options.nc != None):
        configFile = options.nc

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
    shutil.copyfile(dep["viewer.py"], os.path.join(
        io["outputDir"], "viewer.py"))
    shutil.copyfile(dep["plots.py"], os.path.join(io["outputDir"], "plots.py"))

    return dep, io, opt, var, prop, inte


def plyRun(dep, io, opt, var, prop, inte):
    '''function run the driver function starting with .ply mesh'''
    # create starting mesh
    if io["generateGeometry"] == True:
        pymem3dg.genIcosphere(nSub=io["nSub"], path=io["refMesh"], R=io["R"])
    # elif io["inputMesh"] == "UVsphere.ply":
    #     pymem3dg.genUVsphere( nSub = io["nSub"]
    #                     , path = io["inputMesh"])

    # run simulation
    return pymem3dg.driver_ply(verbosity=io["verbosity"],
                               inputMesh=io["inputMesh"],
                               outputDir=io["outputDir"],
                               refMesh=io["refMesh"],
                               nSub=io["nSub"],

                               isVertexShift=opt["isVertexShift"],
                               isProtein=opt["isProtein"],
                               isLocalCurvature=opt["isLocalCurvature"],
                               isReducedVolume=opt["isReducedVolume"],

                               H0=var["H0*R"] / io["R"],
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
                               integration=inte["method"],
                               isBacktrack=inte["options"]["isBacktrack"],
                               rho=inte["options"]["rho"],
                               c1=inte["options"]["c1"],
                               ctol=inte["options"]["ctol"],
                               isAugmentedLagrangian=inte["options"]["isAugmentedLagrangian"])


def ncRun(dep, io, opt, var, prop, inte):
    '''function run the driver function starting with data stored in netcdf .nc file'''
    return pymem3dg.driver_nc(verbosity=io["verbosity"],
                              trajFile=io["trajFile"],
                              startingFrame=io["startingFrame"],
                              outputDir=io["outputDir"],

                              isVertexShift=opt["isVertexShift"],
                              isProtein=opt["isProtein"],
                              isLocalCurvature=opt["isLocalCurvature"],
                              isReducedVolume=opt["isReducedVolume"],

                              H0=var["H0*R"] / io["R"],
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
                              integration=inte["method"],
                              isBacktrack=inte["options"]["isBacktrack"],
                              rho=inte["options"]["rho"],
                              c1=inte["options"]["c1"],
                              ctol=inte["options"]["ctol"],
                              isAugmentedLagrangian=inte["options"]["isAugmentedLagrangian"])

def genPlots(io):
    '''function that generate the .png plots for netcdf trajectory file'''
    sys.path.insert(1, io["outputDir"])
    from plots import plot
    plot(trajnc=os.fspath(os.path.abspath(Path(io['outputDir'] + "/traj.nc"))),
         figureFile=os.fspath(os.path.abspath(
             Path(io['outputDir'] + "/traj_plot.pdf"))),
         show=False, save=True)

def plySystem(dep, io, opt, var, prop, inte):
    '''function run the driver function starting with .ply mesh'''
    # create starting mesh
    if io["generateGeometry"] == True:
        pymem3dg.genIcosphere(nSub=io["nSub"], path=io["refMesh"], R=io["R"])
    # elif io["inputMesh"] == "UVsphere.ply":
    #     pymem3dg.genUVsphere( nSub = io["nSub"]
    #                     , path = io["inputMesh"])

    # run simulation
    return pymem3dg.system_ply(verbosity=io["verbosity"],
                               inputMesh=io["inputMesh"],
                               outputDir=io["outputDir"],
                               refMesh=io["refMesh"],
                               nSub=io["nSub"],

                               isVertexShift=opt["isVertexShift"],
                               isProtein=opt["isProtein"],
                               isLocalCurvature=opt["isLocalCurvature"],
                               isReducedVolume=opt["isReducedVolume"],

                               H0=var["H0*R"] / io["R"],
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
                               integration=inte["method"],
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
        plyRun(dep, io, opt, var, prop, inte)
    elif(options.nc != None):
        ncRun(dep, io, opt, var, prop, inte)

    # testing python binding 
    # f = plySystem(dep, io, opt, var, prop, inte)
    # print(f)
    # print(f.insidePressure)
    # # import pdb
    # # pdb.set_trace()
    # f.getBindingForces()

    # generate plots based on netcdf trajectory
    genPlots(io)
