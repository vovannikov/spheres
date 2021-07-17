import os
import numpy as np
import multiprocessing
from os import listdir
from os.path import isfile, join
import csv
import multiprocessing
import time
import argparse
from pathlib import Path

# import settings
from settings import *

# mesh template
meshTpl = "universal"

# Do real simulations of just debug
doRealSimulations = True

def runCaseWrapper(data):
    comb = data['comb']
    options = data['options']
    pathMsh = data['pathMsh']
    pathLog = data['pathLog']
    pathPrm = data['pathPrm']

    runCase(comb, options, pathMsh, pathLog, pathPrm)

def runCase(comb, options, pathMsh, pathLog, pathPrm):

    geoSym = "3" if options['geometry'] == "3d" else "2"

    procId = multiprocessing.current_process()
    print("#{}: dA={}, dB={}, nratio={}, E={}".format(procId, comb['dA'], comb['dB'], comb['rNrat'], comb['E']))

    R1 = comb['dA'] / 2.
    R2 = comb['dB'] / 2.
    Rn = min(R1, R2) * comb['rNrat']
    E = comb['E'] * 1e-3

    h1 = np.sqrt(R1*R1 - Rn*Rn)
    h2 = np.sqrt(R2*R2 - Rn*Rn)
    offset = h1 + h2

    strH1 = "{:.10f}".format(h1)
    strH2 = "{:.10f}".format(h2)

    caseFileName = "spheres{}_R1={:.1f}_R2={:.1f}_Rn={:.4f}_E={:.0f}".format(geoSym, R1, R2, Rn, comb['E'])
    fileMsh = "{}/{}.msh".format(pathMsh, caseFileName)
    fileLog = "{}/{}.log".format(pathLog, caseFileName)
    fileSet = "{}/{}.prm".format(pathPrm, caseFileName)

    # If log file exists, then we skip this case
    if os.path.isfile(fileLog) and not(options['force']):
        return

    # Generate mesh with gmsh
    pathGeo = os.path.join(pathSpheresSrc, "settings/spheres{}-{}-tpl.geo".format(geoSym, meshTpl))
    if doRealSimulations and (not(os.path.isfile(fileMsh)) or options['mesh']):
        cmdGmsh = "gmsh -{} {} -setnumber R1 {} -setnumber R2 {} -setnumber Rn {} -setnumber h1 {} -setnumber h2 {} -o {}".format(geoSym, pathGeo, R1, R2, Rn, strH1, strH2, fileMsh)
        stream = os.popen(cmdGmsh)
        output = stream.read()
        #os.system(cmdGmsh)

    offset = float(strH1) + float(strH2)

    fileParamsTpl = open(pathParams, "r")
    dataParamsTpl = fileParamsTpl.read()

    caseParams = dataParamsTpl
    caseParams = caseParams.replace("%mesh%", fileMsh)
    caseParams = caseParams.replace("%stiffness%", options['stiffness'])
    caseParams = caseParams.replace("%x2%", str(offset))
    caseParams = caseParams.replace("%E%", str(E))
    caseParams = caseParams.replace("%log%", fileLog)
    caseParams = caseParams.replace("%right%", str(options['right']).lower())

    with open(fileSet, "w") as fileCaseSettings:
        fileCaseSettings.write(caseParams)

    # Run spheres executable
    if doRealSimulations:
        cmdSpheres = "{} {} {}d".format(pathExec, fileSet, geoSym)
        
        # This does not work properly in MPI
        stream = os.popen(cmdSpheres)
        output = stream.read()
        #os.system(cmdSpheres)
        #icomm = comm.Spawn(pathExec, args=fileSet, maxprocs=1, root=0)

def runAllCases(options):

    # 2D vs 3D simulations:
    geoSym = "3" if options['geometry'] == "3d" else "2"
    secSym = "right" if options['right'] else "left"

    # What stiffness do we analyze
    #stiffnessMode = "tension"
    #stiffnessMode = "torsion"
    #stiffnessMode = "bending_displacement"
    #stiffnessMode = "bending_rotation"
    stiffnessMode = options['stiffness']

    # General folders
    pathCase = os.path.join(pathJobs, "{}/{}d/{}".format(stiffnessMode, geoSym, secSym))

    if not(meshTpl == stiffnessMode):
        pathMsh = os.path.join(pathBatch, "mesh")
    else:
        pathMsh = os.path.join(pathCase, "mesh")

    pathLog = os.path.join(pathCase, "log")
    pathPrm = os.path.join(pathCase, "prm")

    # Check and create directories
    Path(pathCase).mkdir(parents=True, exist_ok=True)
    Path(pathMsh).mkdir(parents=True, exist_ok=True)
    Path(pathLog).mkdir(parents=True, exist_ok=True)
    Path(pathPrm).mkdir(parents=True, exist_ok=True)

    nproc = 6

    time_start = time.time()

    # List of params to vary
    lstD = [38, 45, 63, 75, 90, 106, 125, 150, 175, 200] # mkm, total = 10
    lstRnRatio = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5] # mkm, total = 9
    lstE = [5, 10, 20, 40, 70, 100, 130, 160, 200] # mkm, total = 9

    # Debug reduced list
    lstD = [38, 45]
    lstRnRatio = [0.01, 0.02]
    lstE = [5, 10]

    nDiam = int(len(lstD) * (len(lstD) + 1) / 2)
    nRatio = len(lstRnRatio)
    nE = len(lstE)

    nTotal = nDiam * nRatio * nE

    print("# of diameters combinations = {}".format(nDiam))
    print("# of neck ratios = {}".format(nRatio))
    print("# of Young's modulii = {}".format(nE))
    print("# of combinations = {}".format(nTotal))

    lstCombinations = []

    for dA in lstD:
        for dB in lstD:
            if dA >= dB:
                for rNratio in lstRnRatio:
                    for E in lstE:
                        lstCombinations.append({
                            'comb': {'dA': dA, 'dB': dB, 'rNrat': rNratio, 'E': E},
                            'options': options,
                            'pathMsh': pathMsh, 
                            'pathLog': pathLog, 
                            'pathPrm': pathPrm
                        })

    with multiprocessing.Pool(processes=nproc) as pool:
        pool.map(runCaseWrapper, lstCombinations)

    time_end = time.time()

    hours, rem = divmod(time_end-time_start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("Wall time: {:0>2}h:{:0>2}m:{:05.2f}s".format(int(hours),int(minutes),seconds))

    # Now merge all logs int a single data file
    logFiles = [f for f in listdir(pathLog) if isfile(join(pathLog, f))]
    logFiles.sort()

    if stiffnessMode == "tension":
        extraFields = ["Kuz"]
    elif stiffnessMode == "torsion":
        extraFields = ["Krz"]
    elif stiffnessMode == "bending_displacement":
        if options['right']:
            extraFields = ["Kb13", "Kb23", "Kb33", "Kb43"]
        else:
            extraFields = ["Kb11", "Kb21", "Kb31", "Kb41"]
    elif stiffnessMode == "bending_rotation":
        if options['right']:
            extraFields = ["Kb14", "Kb24", "Kb34", "Kb44"]
        else:
            extraFields = ["Kb12", "Kb22", "Kb32", "Kb42"]

    csvMain = ["R1", "R2", "Rn", "E"]
    csvHeader = csvMain + extraFields
    lstRows = []

    for fileCase in logFiles:
        fileCaseFull = os.path.join(pathLog, fileCase)

        fileCaseNoExt = os.path.splitext(fileCase)[0]
        fparts = fileCaseNoExt.split('_')
        csvRow = [0] * (len(csvMain))

        for part in fparts:
            pair = part.split('=')
            if len(pair) == 2:
                pname = pair[0]
                pvalue = pair[1]
        
                pid = csvMain.index(pname)
                csvRow[pid] = pvalue

        with open(fileCaseFull, 'r') as f:
            reader = csv.reader(f, delimiter=' ', skipinitialspace=True)

            loadMagnitude = 1.0
            csvStiffness = [0] * (len(extraFields))
            for row in reader:
                if row[0] == "time":
                    continue
                
                t = float(row[0])
                if t > 0.9:
                    loadMagnitude = float(row[1])
                    for i in range(len(extraFields)):
                        try:
                            csvStiffness[i] = float(row[2 + i]) / loadMagnitude
                        except ValueError:
                            csvStiffness[i] = 0

            csvRow += csvStiffness

        lstRows.append(csvRow)

    # Write resultant CSV file
    csvStifFile = os.path.join(pathCase, "{}_stiffness_{}d.csv".format(stiffnessMode, geoSym))
    with open(csvStifFile, 'w') as myfile:
        wr = csv.writer(myfile, delimiter=',')
        wr.writerow(csvHeader)

        [wr.writerow(row) for row in lstRows]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Simulation parameters')
    parser.add_argument("-g", "--geometry", help="Geometry type", choices=['2d', '3d'], required=True)
    parser.add_argument("-s", "--stiffness", help="Case to analyze", required=True,
        choices=['tension', 'torsion', 'bending_displacement', 'bending_rotation'])
    parser.add_argument("-m", "--mesh", help="Force mesh update", action='store_true', default=False)
    parser.add_argument("-f", "--force", help="Force computations", action='store_true', default=False)
    parser.add_argument("-r", "--right", help="Is right section active", action='store_true', default=True)

    args = parser.parse_args()

    options = {
        'geometry': args.geometry,
        'stiffness': args.stiffness,
        'mesh': args.mesh,
        'force': args.force,
        'right': args.right,
    }

    runAllCases(options)

