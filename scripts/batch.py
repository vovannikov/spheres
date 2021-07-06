import os
import numpy as np
from mpi4py import MPI
from os import listdir
from os.path import isfile, join
import csv

# Do real simulations of just debug
doRealSimulations = True

# paths
pathGeo = "/home/kapusta/development/spheres/settings/spheres2-tpl.geo"
pathParams = "/home/kapusta/development/spheres/settings/settings-tpl.prm"
pathBatch = "/home/kapusta/work/batch"
pathExec = "/home/kapusta/development/spheres-release/src/spheres_example"

pathMsh = pathBatch + "/mesh"
pathLog = pathBatch + "/log"
pathSet = pathBatch + "/prm"

# MPI data
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def runCase(icase, comb):

    rank = comm.Get_rank()
    print("# rank={} case={}: dA={}, dB={}, nratio={}, E={}".format(rank, icase, comb['dA'], comb['dB'], comb['rNrat'], comb['E']))

    R1 = comb['dA'] / 2.
    R2 = comb['dB'] / 2.
    Rn = min(R1, R2) * comb['rNrat']
    E = comb['E'] * 1e-3

    h1 = np.sqrt(R1*R1 - Rn*Rn)
    h2 = np.sqrt(R2*R2 - Rn*Rn)
    offset = h1 + h2

    strH1 = "{:.10f}".format(h1)
    strH2 = "{:.10f}".format(h2)

    caseFileName = "spheres2_R1={}_R2={}_Rn={}_E={}".format(R1, R2, Rn, E)
    fileMsh = "{}/{}.msh".format(pathMsh, caseFileName)
    fileLog = "{}/{}.log".format(pathLog, caseFileName)
    fileSet = "{}/{}.prm".format(pathSet, caseFileName)

    # Generate mesh with gmsh
    if doRealSimulations:
        cmdGmsh = "gmsh -2 {} -setnumber R1 {} -setnumber R2 {} -setnumber Rn {} -setnumber h1 {} -setnumber h2 {} -o {}".format(pathGeo, R1, R2, Rn, strH1, strH2, fileMsh)
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
        cmdSpheres = "{} {}".format(pathExec, fileSet)
        
        # This does not work properly in MPI
        stream = os.popen(cmdSpheres)
        output = stream.read()
        #os.system(cmdSpheres)
        #icomm = comm.Spawn(pathExec, args=fileSet, maxprocs=1, root=0)

# This is to be available on all processors

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

if rank == 0:

    print("# of diameters combinations = {}".format(nDiam))
    print("# of neck ratios = {}".format(nRatio))
    print("# of Young's modulii = {}".format(nE))
    print("# of combinations = {}".format(nTotal))

    lstCombinations = []

    for dA in lstD:
        for dB in lstD:
            if dB >= dA:
                for rNratio in lstRnRatio:
                    for E in lstE:
                        lstCombinations.append({'dA': dA, 'dB': dB, 'rNrat': rNratio, 'E': E})

    icase = 0
    if size == 1:
        for comb in lstCombinations:
            runCase(icase, comb)
            icase += 1
    else:
        for comb in lstCombinations:

            nextReady = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)

            data = {'icase':icase, 'comb': comb}
            req = comm.send(data, dest=nextReady, tag=0)
            icase += 1

        for proc in range(1, size):
            data = {'icase':nTotal, 'comb': None}
            comm.send(data, dest=proc, tag=0)

else:
    
    while (True):

        comm.send(rank, dest=0, tag=1)

        data = comm.recv(source=0)
        icase = data['icase']
        if icase < nTotal:
            comb = data['comb']
            runCase(icase, comb)
        else:
            print("Proc {} exiting ...".format(rank))
            break

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
                stiffness = float(row[1])
            except ValueError:
                stiffness = 0

        csvRow[-1] = stiffness

    lstRows.append(csvRow)

# Write resultant CSV file
csvStifFile = os.path.join(pathBatch, "tensile_stiffness_2d_coarse.csv")
with open(csvStifFile, 'w') as myfile:
    wr = csv.writer(myfile, delimiter=',')
    wr.writerow(csvHeader)

    [wr.writerow(row) for row in lstRows]

