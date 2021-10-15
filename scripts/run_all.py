import os
import shutil
import numpy as np
from batch import runAllCases
from settings import *

if doCleanup:
    shutil.rmtree(pathJobs)
    shutil.rmtree(pathMesh)

lstJobs = [
    {'stiffness': 'tension', 'right': True},
    {'stiffness': 'torsion', 'right': True},
    {'stiffness': 'bending_displacement', 'right': False},
    {'stiffness': 'bending_rotation', 'right': False},
    {'stiffness': 'bending_displacement', 'right': True},
    {'stiffness': 'bending_rotation', 'right': True}
]

for job in lstJobs:

    options = {
        'geometry': geometry,
        'stiffness': job['stiffness'],
        'mesh': doForceMesh,
        'force': doForceCalc,
        'right': job['right'],
    }

    print("Running job geometry={}, stiffness = {}, forceMesh = {}, forceCalc = {}, isRight = {}".format(geometry, job['stiffness'], doForceMesh, doForceCalc, job['right']))

    runAllCases(options)


# Now merge all csv together
lstCsv = []
for root, dirs, files in os.walk(pathJobs):
    for file in files:
        if file.endswith(".csv"):
            if not(os.path.exists(os.path.join(pathJobs, file))):
                filePath = os.path.join(root, file)
                lstCsv.append(filePath)

resPre = []
resData = []
isPreDone = False
for csvFile in lstCsv:

    resTable = np.genfromtxt(csvFile, delimiter=',', names=True)

    for cname in resTable.dtype.names:
        if cname.startswith('K'):
            resData.append({'header': cname, 'data': resTable[cname]})
        elif not(isPreDone):
            resPre.append({'header': cname, 'data': resTable[cname]})

    isPreDone = True

def mycmp(item1, item2):
    if item1['header'] < item2['header']:
        return -1
    elif item1['header'] > item2['header']:
        return 1
    else:
        return 0

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0  
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

resData.sort(key=cmp_to_key(mycmp))

header = []
csv = []

for col in resPre:
    header.append(col['header'])
    csv.append(col['data'])

for col in resData:
    header.append(col['header'])
    csv.append(col['data'])

# Resultant file name where everything is merged together
pathCsvResultant = os.path.join(pathJobs, resultantFileName)

if os.path.exists(pathCsvResultant):
    os.remove(pathCsvResultant)

if len(csv[0].shape) == 0:
    csv = np.array([c.item() for c in csv])
    csv = csv.reshape(1, csv.shape[0])
else:
    csv = np.transpose(csv)

np.savetxt(pathCsvResultant, csv, header=','.join(header), comments='', delimiter=',')
