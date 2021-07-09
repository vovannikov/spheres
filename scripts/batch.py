import os
import numpy as np
import multiprocessing
from os import listdir
from os.path import isfile, join
import csv
import multiprocessing
import time

# 2D vs 3D simulations:
is3D = True
geoSym = "3" if is3D else "2"

# Do real simulations of just debug
doRealSimulations = True

# paths
pathParams = "/home/kapusta/development/spheres/settings/settings-tpl.prm"
pathExec = "/home/kapusta/development/spheres-release/src/spheres_runner"
pathGeo = "/home/kapusta/development/spheres/settings/spheres{}-tpl.geo".format(geoSym)
pathBatch = "/home/kapusta/work/batch/{}d".format(geoSym)

pathMsh = pathBatch + "/mesh"
pathLog = pathBatch + "/log"
pathSet = pathBatch + "/prm"

def runCase(comb):

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
    fileSet = "{}/{}.prm".format(pathSet, caseFileName)

    # If log file exists, then we skip this case
    if os.path.isfile(fileLog):
        return

    # Generate mesh with gmsh
    if doRealSimulations:
        cmdGmsh = "gmsh -{} {} -setnumber R1 {} -setnumber R2 {} -setnumber Rn {} -setnumber h1 {} -setnumber h2 {} -o {}".format(geoSym, pathGeo, R1, R2, Rn, strH1, strH2, fileMsh)
        stream = os.popen(cmdGmsh)
        output = stream.read()
        #os.system(cmdGmsh)

    offset = float(strH1) + float(strH2)

    fileParamsTpl = open(pathParams, "r")
    dataParamsTpl = fileParamsTpl.read()

    caseParams = dataParamsTpl
    caseParams = caseParams.replace("%mesh%", fileMsh)
    caseParams = caseParams.replace("%x2%", str(offset))
    caseParams = caseParams.replace("%E%", str(E))
    caseParams = caseParams.replace("%log%", fileLog)

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

nproc = 6

if __name__ == '__main__':

    time_start = time.time()

    # List of params to vary
    lstD = [38, 45, 63, 75, 90, 106, 125, 150, 175, 200] # mkm, total = 10
    lstRnRatio = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5] # mkm, total = 9
    lstE = [5, 10, 20, 40, 70, 100, 130, 160, 200] # mkm, total = 9

    # Debug reduced list
    #lstD = [38, 45]
    #lstRnRatio = [0.01, 0.02]
    #lstE = [5, 10]

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
                        lstCombinations.append({'dA': dA, 'dB': dB, 'rNrat': rNratio, 'E': E})

    with multiprocessing.Pool(processes=nproc) as pool:
        pool.map(runCase, lstCombinations)

    time_end = time.time()

    hours, rem = divmod(time_end-time_start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("Wall time: {:0>2}h:{:0>2}m:{:05.2f}s".format(int(hours),int(minutes),seconds))

# Now merge all logs int a single data file
logFiles = [f for f in listdir(pathLog) if isfile(join(pathLog, f))]

csvHeader = ["R1", "R2", "Rn", "E", "Kzz"]
lstRows = []

for fileCase in logFiles:
    fileCaseFull = os.path.join(pathLog, fileCase)

    fileCaseNoExt = os.path.splitext(fileCase)[0]
    fparts = fileCaseNoExt.split('_')
    csvRow = [0] * (len(csvHeader))

    for part in fparts:
        pair = part.split('=')
        if len(pair) == 2:
            pname = pair[0]
            pvalue = pair[1]
    
            pid = csvHeader.index(pname)
            csvRow[pid] = pvalue

    with open(fileCaseFull, 'r') as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)

        stiffness = 0
        for row in reader:
            try:
                t = float(row[0])
                if t > 0.9:
                    stiffness = float(row[1])
            except ValueError:
                stiffness = 0

        csvRow[-1] = stiffness

    lstRows.append(csvRow)

# Write resultant CSV file
csvStifFile = os.path.join(pathBatch, "tensile_stiffness_{}d.csv".format(geoSym))
with open(csvStifFile, 'w') as myfile:
    wr = csv.writer(myfile, delimiter=',')
    wr.writerow(csvHeader)

    [wr.writerow(row) for row in lstRows]

