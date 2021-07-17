import os

# paths TO DEFINE
pathHome = os.path.expanduser("~")
pathSpheresSrc = os.path.join(pathHome, "development/spheres")
pathExec = os.path.join(pathHome, "development/spheres-release/src/spheres_runner")
pathBatch = os.path.join(pathHome, "work/spheres")

resultantFileName = 'stiffness_3d.csv'

# Dependent paths, do not touch this
pathParams = os.path.join(pathSpheresSrc, "settings/settings-tpl.prm")
pathJobs = os.path.join(pathBatch, "jobs")