#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from glob import glob
from sys import exit


wd = "/home/saultyevil/Dropbox/DiskWinds/PySims/tde/sanity_checks/cv_cno_f1e-2"
wd = "./"
root = "tde_cv"

glob_directory = "{}/diag_{}/{}_*.diag".format(wd, root, root)
diag_files = glob(glob_directory)
diag_files = sorted(diag_files)
# print(diag_files)

broken_diag = []

ndiag = len(diag_files)

l_adiabatic = []
l_low_freq = []
l_absorbed = []
l_scatters = []
l_star = []
l_disk = []
l_errors = []
l_unknown = []

for i in range(ndiag):

    diag = diag_files[i]

    try:
        with open(diag, "r") as f:
            lines = f.readlines()
    except IOError:
        print("couldn't open diag file {}".format(diag))
        broken_diag.append(diag)
        continue

    for line in lines:

        j = 0

        # TODO: this can be optimised AHHHHHHHH

        # print(line.split())

        # print(i)

        if i == 0:
            if line.find("!!python: luminosity lost by adiabatic kpkt destruction") != -1:
                l_adiabatic.append(float(line.split()[7]))
            if line.find("!!python: luminosity lost to low-frequency free-free") != -1:
                l_low_freq.append(float(line.split()[6]))
            if line.find("!!python: luminosity lost by being completely absorbed") != -1:
                l_absorbed.append(float(line.split()[7]))
            if line.find("!!python: luminosity lost by too many scatters") != -1:
                l_scatters.append(float(line.split()[7]))
            if line.find("!!python: luminosity lost by hitting the star") != -1:
                l_star.append(float(line.split()[7]))
            if line.find("!!python: luminosity lost by hitting the disk") != -1:
                l_disk.append(float(line.split()[7]))
            if line.find("!!python: luminosity lost by errors") != -1:
                l_errors.append(float(line.split()[5]))
            if line.find("!!python: luminosity lost by the unknown") != -1:
                l_unknown.append(float(line.split()[6]))
        else:
            if line.find("!!python: luminosity lost by adiabatic kpkt destruction") != -1:
                l_adiabatic[j] += (float(line.split()[7]))
            if line.find("!!python: luminosity lost to low-frequency free-free") != -1:
                l_low_freq[j] += (float(line.split()[6]))
            if line.find("!!python: luminosity lost by being completely absorbed") != -1:
                l_absorbed[j] += (float(line.split()[7]))
            if line.find("!!python: luminosity lost by too many scatters") != -1:
                l_scatters[j] += (float(line.split()[7]))
            if line.find("!!python: luminosity lost by hitting the star") != -1:
                l_star[j] += (float(line.split()[7]))
            if line.find("!!python: luminosity lost by hitting the disk") != -1:
                l_disk[j] += (float(line.split()[7]))
            if line.find("!!python: luminosity lost by errors") != -1:
                l_errors[j] += (float(line.split()[5]))
            if line.find("!!python: luminosity lost by the unknown") != -1:
                l_unknown[j] += (float(line.split()[6]))
                j += 1

print("!!python: luminosity lost by adiabatic kpkt destruction", l_adiabatic)
print("!!python: luminosity lost to low-frequency free-free", l_low_freq)
print("!!python: luminosity lost by being completely absorbed", l_absorbed)
print("!!python: luminosity lost by too many scatters", l_scatters)
print("!!python: luminosity lost by hitting the star", l_star)
print("!!python: luminosity lost by hitting the disk", l_disk)
print("!!python: luminosity lost by errors", l_errors)
print("!!python: luminosity lost by the unknown", l_unknown)