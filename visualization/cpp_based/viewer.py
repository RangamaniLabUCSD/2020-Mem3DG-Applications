import sys
import json
import os
import optparse
import pymem3dg
from pathlib import Path

# parse the config file and options
optparser = optparse.OptionParser()
optparser.add_option(
    "-a",
    "--animate",
    dest="animate",
    type="string",
    help="visualize animation of trajectory file",
)
optparser.add_option(
    "-s",
    "--snapshot",
    dest="snapshot",
    type="string",
    help="visualize snapshot geometric file",
)
(options, args) = optparser.parse_args()

# make video directory if not exist
cwd = os.getcwd()
videoDir = os.path.join(cwd, "./video")
if not os.path.exists(videoDir):
    os.mkdir(videoDir)

####################################################
# initialize option for visualization, CHANGE HERE #
####################################################
# Quantities data needed for visualization
Q = pymem3dg.Quantities()
Q.ref_coord = True
Q.velocity = True
Q.mean_curvature = True
Q.gauss_curvature = True
Q.spon_curvature = True
Q.ext_force = True
Q.physical_force = True
Q.capillary_force = True
Q.osmotic_force = True
Q.bending_force = True
Q.line_force = True
Q.mask = False
Q.H_H0 = True
Q.the_point = False
Q.smoothing_mask = False
Q.chemical_potential = True
Q.bending_potential = True
Q.diffusion_potential = True
Q.adsorption_potential = True

# GUI & misc: optional arguments for viewers
transparency = 1
angle = 0
fov = 50
edgeWidth = 1
isShow = True
isSave = False
screenshotName = "screenshot.png"

######################################################
# run viewer based on option arguments and extension #
######################################################
if options.snapshot != None:
    options.snapshot = os.fspath(os.path.abspath(Path(options.snapshot)))
    ext = os.path.splitext(options.snapshot)[1][1:]
    # run viewer
    if ext == "nc":
        pymem3dg.snapshot_nc(
            fileName=options.snapshot,
            options=Q,
            frame=options[0],
            transparency=transparency,
            angle=angle,
            fov=fov,
            edgeWidth=edgeWidth,
            isShow=isShow,
            isSave=isSave,
            screenshotName=screenshotName,
        )

    elif ext == "ply":
        pymem3dg.snapshot_ply(
            fileName=options.snapshot,
            options=Q,
            transparency=transparency,
            fov=fov,
            edgeWidth=edgeWidth,
        )

elif options.animate != None:
    options.animate = os.fspath(os.path.abspath(Path(options.animate)))
    ext = os.path.splitext(options.animate)[1][1:]
    # run viewer
    if ext == "nc":
        pymem3dg.animate_nc(
            fileName=options.animate,
            options=Q,
            transparency=transparency,
            fov=fov,
            edgeWidth=edgeWidth,
        )
    else:
        pymem3dg.animate_ply(
            framesDir=options.animate,
            options=Q,
            frameNum=[int(args[0]), int(args[1])],
            transparency=transparency,
            fov=fov,
            edgeWidth=edgeWidth,
        )
