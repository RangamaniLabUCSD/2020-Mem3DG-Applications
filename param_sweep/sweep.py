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

                H0 = 0,
                Vt = 0.65,
                ptInd = 0,  
                Kf = 0,
                conc = 25,
                height = 0,
                
                Kb  = kb,
                Kse = 0,     
                Ksl = 0.1,
                Kst = 0.01,		
                Ksg = [15,100],		
                Kv  = [15,100], 	
                kt = 0, 
                gamma = 10, 
                
                h = 0.0001,
                T = 20,
                eps = 0.02,
                closeZone = 10000,
                increment = 0.1,
                tSave = 10,
                tMollify = 10)
    
def runSims():
    jobs = []
    for i, kb in enumerate(np.arange(0,0.1,0.01)):
        for replicate in np.arange(0,10000):
            path = f'run_{i}_{replicate}/'
            jobs.append((kb,path,))
    r = process_map(worker, jobs, max_workers=4)

if __name__ == "__main__":
    runSims()
