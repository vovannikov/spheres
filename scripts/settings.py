import os

# paths TO DEFINE
pathHome = os.path.expanduser("~")
pathSpheresSrc = os.path.join(pathHome, "development/spheres")
pathExec = os.path.join(pathHome, "development/spheres-release/src/spheres_runner")
pathBatch = os.path.join(pathHome, "work/spheres")

# Number of processes
nproc = 6

# Remote 
resultantFileName = 'stiffness_3d.csv'

# Dependent paths, do not touch this !!!!
pathParams = os.path.join(pathSpheresSrc, "settings/settings-tpl.prm")
pathJobs = os.path.join(pathBatch, "jobs")
pathMesh = os.path.join(pathBatch, "mesh")

# Settings for ran_all script
geometry = '3d'
doForceMesh = False
doForceCalc = False
doCleanup = True
doSaveVtk = True

# Mesh quality
meshLc = 0.5 # 2.0 - super coarse for debug, 0.5 - coarse computations, 0.1 - ok computations, 0.02 - fine computations

# =======================
# Lists of params to vary

# Complete list
lstD = [38, 45, 63, 75, 90, 106, 125, 150, 175, 200, 250, 300] # mkm, total = 12
lstRnRatio = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6] # mkm, total = 10
lstE = [5, 10, 25, 50, 100, 150, 200] # mkm, total = 7

# Debug reduced list
#lstD = [38, 45]
##lstRnRatio = [0.01, 0.02]
##lstE = [5, 10]
#lstRnRatio = [0.4]
#lstE = [5]