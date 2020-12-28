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
      pymem3dg.view_animation(options.shape,ref_coord = 0, velocity = 0,
                          mean_curvature = 0,  spon_curvature = 0,
                        ext_pressure = 0, physical_pressure = 0,
                        capillary_pressure = 0,
                        bending_pressure = 0, line_pressure = 0, mask = 0, H_H0 = 0)
  elif ext == "ply": 
      pymem3dg.viewer(options.shape)
elif(options.all != None):
  ext = os.path.splitext(options.all)[1][1:]
  # run viewer 
  if ext == "nc":
      pymem3dg.view_animation(options.all, ref_coord = 1, velocity = 1,
                          mean_curvature = 1,  spon_curvature = 1,
                        ext_pressure = 1, physical_pressure = 1,
                        capillary_pressure = 1,
                        bending_pressure = 1, line_pressure = 1, mask = 1, H_H0 = 1)
  elif ext == "ply": 
      pymem3dg.viewer(options.all, mean_curvature = 1, spon_curvature = 1,
          ext_pressure = 1, physical_pressure = 1, capillary_pressure = 1,
            bending_pressure = 1, line_pressure = 1)
