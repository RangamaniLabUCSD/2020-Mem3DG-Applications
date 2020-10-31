#!/bin/python

import sys
import os
from pathlib import Path
import argparse
import shutil
import pymem3dg

import multiprocessing
import subprocess
from tqdm.contrib.concurrent import process_map
import numpy as np

import concurrent.futures

inPath = os.fspath(Path('../input-file/patch.ply'))
refPath = os.fspath(Path('../input-file/patch.ply'))

def worker(args):
    H0 = args[0]
    oPath = args[1]
    # print(oPath)
    if not os.path.exists(oPath):
        os.mkdir(oPath)

    # run simulation
    pymem3dg.driver_ply(inputMesh = inPath,
                outputDir = oPath,
                refMesh = refPath,
                radius = 0.9,
                nSub = 0,

                isTuftedLaplacian = False,
                mollifyFactor = 1e-3,
                isVertexShift = False,
                isProtein = False,

                epsilon = 15e-5,
                Bc = 40,
                H0 = H0,
                sharpness = 100,
                r_H0 = 0.1,
                Vt = 1,
                ptInd = 672,  
                Kf = 0,
                conc = 25,
                height = 0,
                
                Kb  = 8.22e-5,
                Kse = 0,     
                Ksl = 0.05,
                Kst = 0,		
                Ksg = [0.1,0.1],		
                Kv  = [5e-2,5e-2], 	
                kt = 0, 
                gamma = 2, 
                
                h = 1e-5,
                T = 50,
                eps = 0.002,
                closeZone = 1000,
                increment = 0.1,
                tSave = 1e-1,
                tMollify = 100)
    
def runSims():
    jobs = []
    for i, H0 in enumerate(np.arange(35,45,1)):
        for replicate in np.arange(0,1):
            path = f'run_{i}_{replicate}/'
            jobs.append((H0,path))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(worker, jobs)     

    # r = process_map(worker, jobs, max_workers=4)

if __name__ == "__main__":
    runSims()
