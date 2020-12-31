import sys
import json
import os
import optparse
import pymem3dg
# parse the config file and options
optparser = optparse.OptionParser()
optparser.add_option('-s', '--shape', dest="shape", type="string",
                     help="visualization file in .nc/.ply format")
optparser.add_option('-a', '--all', dest="all", type="string",
                     help="visualization file in .nc/.ply format")
(options, args) = optparser.parse_args()
# argparser = argparse.ArgumentParser()
# parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
# args = argparser.parse_args()

# make video directory if not exist
cwd = os.getcwd()
videoDir = os.path.join(cwd, './video')
if not os.path.exists(videoDir):
    os.mkdir(videoDir)

if (options.shape != None):
    ext = os.path.splitext(options.shape)[1][1:]
    # run viewer
    if ext == "nc":
        pymem3dg.animation_nc(fileName=options.shape, ref_coord=False, velocity=False,
                              mean_curvature=False,  spon_curvature=False,
                              ext_pressure=False, physical_pressure=False,
                              capillary_pressure=False,
                              bending_pressure=False, line_pressure=False, mask=False, H_H0=False)
    elif ext == "ply":
        pymem3dg.viewer_ply(fileName=options.shape, mean_curvature=False, spon_curvature=False,
                            ext_pressure=False, physical_pressure=False, capillary_pressure=False,
                            bending_pressure=False, line_pressure=False)

elif(options.all != None):
    ext = os.path.splitext(options.all)[1][1:]
    # run viewer
    if ext == "nc":
        pymem3dg.animation_nc(fileName=options.all, ref_coord=True, velocity=True,
                              mean_curvature=True,  spon_curvature=True,
                              ext_pressure=True, physical_pressure=True,
                              capillary_pressure=True,
                              bending_pressure=True, line_pressure=True, mask=True, H_H0=True)

    elif ext == "ply":
        pymem3dg.viewer_ply(fileName=options.all, mean_curvature=True, spon_curvature=True,
                            ext_pressure=True, physical_pressure=True, capillary_pressure=True,
                            bending_pressure=True, line_pressure=True)
