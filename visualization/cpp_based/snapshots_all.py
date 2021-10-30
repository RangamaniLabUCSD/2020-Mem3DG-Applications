import os
import glob
import pymem3dg
import netCDF4 as nc

if __name__ == "__main__":
    # get the working directory
    cwd = os.getcwd()
    cwd_string = os.fspath(os.path.abspath(cwd))

    # recursively parse every trajectory netcdf files
    for filename in glob.iglob(cwd_string + '/**/**.nc', recursive=True):

        ds = nc.Dataset(filename)
        lastFrameIndex = len(ds.dimensions['frame']) - 1

        # extract netcdf file name without extension
        directory, ncFile = os.path.split(filename)
        _, folderName = os.path.split(directory)
        ncFileName = folderName + "_" + os.path.splitext(ncFile)[0]

        pymem3dg.snapshot_nc(fileName=filename, frame=lastFrameIndex, transparency=0.3,
                             angle=0, fov=70, edgeWidth=2, isShow=False, isSave=True,
                             screenshotName=ncFileName + ".png", ref_coord=False, velocity=False,
                             mean_curvature=True,  spon_curvature=False,
                             ext_pressure=False, physical_pressure=False,
                             capillary_pressure=False, inside_pressure=False,
                             bending_pressure=False, line_pressure=False, mask=False, H_H0=False)