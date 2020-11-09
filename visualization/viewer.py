import sys
import json
import argparse
import os

# parse the config file 
parser = argparse.ArgumentParser()
# parser.add_argument("config", help = "configuration file (.json) used for the simulation", type = str)
parser.add_argument("vis_file", help = "Mesh file (.ply, .obj, .etc.) used for visualization, OR \
                    trajectory file (.nc) used for animation", type = str)
args = parser.parse_args()
ext = os.path.splitext(args.vis_file)[1][1:]
# # read the config file
# with open(args.config) as f:
#     data = json.load(f)
# geo = data["parameters"]["geometry"]
# prop = data["parameters"]["properties"]
# var = data["parameters"]["variables"]
# inte = data["parameters"]["integration"]
# dep = data["parameters"]["dependencies"]
# opt = data["parameters"]["options"]

# # find dependency
# sys.path.append(dep["pyddg"])
import pymem3dg

# run viewer 
if ext == "nc":
    pymem3dg.view_animation(args.vis_file, ref_coord = 1, velocity = 1,
                        mean_curvature = 1,  spon_curvature = 1,
                       ext_pressure = 1, physical_pressure = 1,
                       capillary_pressure = 1,
                       bending_pressure = 1)
elif ext == "ply": 
    pymem3dg.viewer(args.vis_file, mean_curvature = 1, spon_curvature = 1,
         ext_pressure = 1, physical_pressure = 1, capillary_pressure = 1,
           bending_pressure = 1)