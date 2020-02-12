#!/usr/bin/python

# Scane a LAMMPS output file and compute the enthalpy of each species

import sys, os
import math
import operator as o
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from scipy import linalg as LA

from dump import dump


def readTypeInfo(lab):
    sour = open(lab + '.lmpdat', 'r')
    for line in sour:
        if line == 'Masses\n':
            break
    mass = []
    sour.next()
    massLine = sour.next()
    while massLine != '\n':
        words = massLine.split()
        mass.append(words[0:2])
        mass[-1][0] = str(mass[-1][0])
        massLine = sour.next()
    for line in sour:
        if line == 'Atoms\n':
            break
    type = []
    sour.next()
    typeLine = sour.next()
    while typeLine != '\n':
        words = typeLine.split()
        type.append([words[0], words[2]])
        typeLine = sour.next()
    sour.close()
    atomMasses = []
    for e in type:
        for element in mass:
            if e[-1] == element[0]:
                atomMass = float(element[-1])
        atomMasses.append(atomMass)
    return atomMasses


def gaussian(list, point, variance):
    def gaussianFunc(x, point, variance):
        return (1 / (2 * math.pi * variance) ** 0.5) * math.exp(-(x - point) ** 2 / 2 / variance)

    for e in list:
        e[1] += gaussianFunc(e[0], point, variance)
    return list


linker = sys.argv[1]
#linker = 'linker205'
inputfile = linker + "-mass-weighted-force-matrix.d"  # Name of LAMMPS file containing position data
outputfile = linker + "-mass-weighted-force-hessian.d"

frame = 0

# First read one frame and then find the number of species

d = dump(inputfile, 0)  # Create dump class but don't read the file all at once

d.map(1, "id", 2, "fx", 3, "fy", 4, "fz")

hessian = [];

time = d.next()
id, fx, fy, fz = d.vecs(time, "id", "fx", "fy", "fz")
d.tselect.none()
d.delete()

nat = len(id)
nmd = nat * 3

fo = np.array([fx, fy, fz]).T
# print fo

for i in range(nmd):
    # for i in range(1):

    time = d.next()
    id, fx, fy, fz = d.vecs(time, "id", "fx", "fy", "fz")
    fp = np.array([fx, fy, fz]).T
    d.tselect.none()
    d.delete()

    time = d.next()
    id, fx, fy, fz = d.vecs(time, "id", "fx", "fy", "fz")
    fm = np.array([fx, fy, fz]).T
    d.tselect.none()
    d.delete()

    K = (fp - fm) * 0.5
    #    K = conv*(fp-fm)/(2.0*delta*mass)
    hessian.append(np.concatenate(K))
# print 'hesian = ',hessian[-1]

hessianPlus = list(np.array(hessian).copy())
for i in range(len(hessian)):
    for j in range(len(hessian[0])):
        hessianPlus[i][j] = -(hessian[i][j] + hessian[j][i]) / 2.0
vals,vecs = LA.eig(hessianPlus)
e_vals = []
for i in range(len(vals)):
    e_vals.append(abs(np.real(vals[i])))
    e_vals[i] = e_vals[i] ** 0.5 * 10 ** -12 / 2 / math.pi
e_vals.sort()
#print(e_vals)
f.close()

vibrationSpectrum = []
for i in np.arange(0, 150, 0.5):
    vibrationSpectrum.append([i, 0])

compressRange = [0, 150]

for e in e_vals:
    vibrationSpectrum = gaussian(vibrationSpectrum, e, 3)
vibrationSpectrum = np.array(vibrationSpectrum)
vibrationSpectrum[:, 1] = vibrationSpectrum[:, 1] / len(vibrationSpectrum[:, 1])
np.save(os.path.join(output_path,linker + ".npy"), vibrationSpectrum[:, 1])
'''
plt.figure()
plt.title(linker+'-Vibrating-frequency-spectrum')
plt.plot(vibrationSpectrum[:, 0], vibrationSpectrum[:, 1])
plt.xlabel('frequency')
plt.ylabel('intensity')
plt.savefig(linker+'-Vibrating-frequency-spectrum.png')
'''
