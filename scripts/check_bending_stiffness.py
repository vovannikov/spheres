from settings import *
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Print bending stiffness structure')
parser.add_argument("-n", "--number", help="Number of entry to print", type=int, required=True)

args = parser.parse_args()

pathCsvResultant = os.path.join(pathJobs, resultantFileName)

csv = np.genfromtxt(pathCsvResultant, delimiter=',', names=True)

n = args.number

K = np.matrix([
    [csv['Kb11'][n], csv['Kb12'][n], csv['Kb13'][n], csv['Kb14'][n]],
    [csv['Kb21'][n], csv['Kb22'][n], csv['Kb23'][n], csv['Kb24'][n]],
    [csv['Kb31'][n], csv['Kb32'][n], csv['Kb33'][n], csv['Kb34'][n]],
    [csv['Kb41'][n], csv['Kb42'][n], csv['Kb43'][n], csv['Kb44'][n]]])

print("Bending stiffness for entry #{}".format(n))
print(K)