#!/bin/python

import sys
import os

import argparse
import shutil
import pymem3dg

import multiprocessing
import subprocess
from tqdm.contrib.concurrent import process_map
import numpy as np

import concurrent.futures

def worker(args):
    kb = args[0]
    oPath = args[1]
    # print(oPath)
    os.mkdir(oPath)

    # run simulation
    pymem3dg.driver(inputMesh = 'slightlyOblate.ply',
                outputDir = oPath,
                refMesh = 'slightlyOblate.ply',
                radius = 10,

                isTuftedLaplacian = False,
                mollifyFactor = 1e-3,
                isVertexShift = False,
                isProtein = true,

                epsilon = 8e-5,
                Bc = 20,
                H0 = 4,
                sharpness = 10,
                r_H0 = 3,
                Vt = 0.65,
                ptInd = 0,  
                Kf = 0,
                conc = 25,
                height = 0,
                
                Kb  = kb,
                Kse = 0,     
                Ksl = 0.01,
                Kst = 0,		
                Ksg = [0.15,1],		
                Kv  = [0.05,0.1], 	
                kt = 0, 
                gamma = 10, 
                
                h = 1e-3,
                T = 1000,
                eps = 0.002,
                closeZone = 1000,
                increment = 0.1,
                tSave = 10,
                tMollify = 100)
    
def runSims():
    jobs = []
    for i, kb in enumerate(np.arange(1e-5,0.1,1e-4)):
        for replicate in np.arange(0,1):
            path = f'run_{i}_{replicate}/'
            jobs.append((kb,path))

    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     results = executor.map(worker, jobs)     

    r = process_map(worker, jobs, max_workers=4)

if __name__ == "__main__":
    runSims()
