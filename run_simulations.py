#!/usr/bin/env python3

import numpy as np
from subprocess import Popen, PIPE

with open("simulations", "r") as f:
    sims = f.readlines()

for sim in sims:
    cd = "cd {}".format(sims)
    print(cd)
    a = Popen(cd, stdout=PIPE, stderr=PIPE, shell=True)
    setup = "Setup_Py_Dir"
    aa = Popen(setup, stdout=PIPE, stderr=PIPE, shell=True)
    pf = "tde.pf"
    run_python = "mpirun py {}".format(pf)
    aaa = Popen(run_python, stdout=PIPE, stderr=PIPE, shell=True)
    convergence = "convergence.py {}".format(pf)
    aaaa = Popen(convergence, stdout=PIPE, stderr=PIPE, shell=True)
    plot = "spec_plot.py tde all -dist 1.079987153448e+27"
    aaaaa = Popen(plot, stdout=PIPE, stderr=PIPE, shell=True)