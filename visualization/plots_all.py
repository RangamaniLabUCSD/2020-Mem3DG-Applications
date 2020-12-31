import os
import glob
from plots import plot

if __name__ == "__main__":
  
    ## get the working directory
    cwd = os.getcwd()
    cwd_string = os.fspath(os.path.abspath(cwd))

    ## recursively parse every trajectory netcdf files 
    for filename in glob.iglob(cwd_string + '/**/**.nc', recursive=True):

        ## extract the foldername and filename
        directory, _ = os.path.split(filename)
        _, folderName = os.path.split(directory)
    
        # Run plot() with input of filename and output as foldername
        plot(trajnc=filename, figureFile=folderName + ".pdf", show=False, save=True)