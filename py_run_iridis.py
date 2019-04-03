#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import argparse
import py_change_parameter
from subprocess import Popen, PIPE

CONVERGED = \
    """
                                                 _
      ___ ___  _ ____   _____ _ __ __ _  ___  __| |
     / __/ _ \| '_ \ \ / / _ \ '__/ _` |/ _ \/ _` |
    | (_| (_) | | | \ V /  __/ | | (_| |  __/ (_| |
     \___\___/|_| |_|\_/ \___|_|  \__, |\___|\__,_|
                                  |___/
    """

NOT_CONVERGED = \
    """
                 _                                                 _
     _ __   ___ | |_    ___ ___  _ ____   _____ _ __ __ _  ___  __| |
    | '_ \ / _ \| __|  / __/ _ \| '_ \ \ / / _ \ '__/ _` |/ _ \/ _` |
    | | | | (_) | |_  | (_| (_) | | | \ V /  __/ | | (_| |  __/ (_| |
    |_| |_|\___/ \__|  \___\___/|_| |_|\_/ \___|_|  \__, |\___|\__,_|
                                                    |___/
    """

ITS_A_MYSTERY = \
    """
     _ _   _                                   _
    (_) |_( )___    __ _   _ __ ___  _   _ ___| |_ ___ _ __ _   _
    | | __|// __|  / _` | | '_ ` _ \| | | / __| __/ _ \ '__| | | |
    | | |_  \__ \ | (_| | | | | | | | |_| \__ \ ||  __/ |  | |_| |
    |_|\__| |___/  \__,_| |_| |_| |_|\__, |___/\__\___|_|   \__, |
                                     |___/                  |___/
    """

CLIM = 0.90


def check_convergence(root: str):
    """
    Check the convergence of a simulation. Assumes the script is being run in
    the correct directory.

    Parameters
    ----------
    root            str
                    The root name of the Python simulation.

    Returns
    -------
    converged       list[str]
                    A list of the convergence fraction for each cycle
    """

    diag_path = "diag_{}/{}_0.diag".format(root, root)
    try:
        with open(diag_path, "r") as f:
            diag = f.readlines()
    except IOError:
        return -1

    converged = []
    for line in diag:
        if line.find("!!Check_converging") != -1:
            tmp1 = line.split()[2].replace("(", "").replace(")", "")
            converged.append(float(tmp1))

    return converged


def run_python(py: str, flags: str, root: str, ncores: str, ion_or_spec: str, verbose: bool = False):
    """
    Run a Python simulation on the Iridis supercomputer.

    If this function is called with "ion" for ion_or_spec, then the number of spectrum
    cycles is set to 0. This is mostly to avoid creating spectra for models which
    have not converged. If the function is called with "spec" then the parameter
    file is updated with spectrum cycles and the photon number is decreased to
    improve run time of the spectrum cycles.

    Parameters
    ----------
    py              str
                    The name of the Python executable
    flags           str
                    Any flags which Python will be called with
    root            str
                    The root name of the Python simulation
    ncores          str
                    The number of cores to run Python with
    ion_or_spec     str
                    An indicating variable to say if simulating an ionisation
                    or a spectral cycle
    verbose         bool, optional
                    Enable verbose logging

    Returns
    -------
    Depending on the convergence of a simulation, a different int will be returned.
    If 0 is returned, then the simulation has reached the desired convergence
    which is set by the variable CLIM. Otherwise, a non-zero int is returned.

    For spec cycles, the function will always return 0.

    """

    print("Running {} cycles for root {}".format(ion_or_spec, root))
    print("CLIM {}".format(CLIM))

    if ion_or_spec == "ion":
        py_change_parameter.change_python_parameter(root, "Spectrum_cycles", "0")
    elif ion_or_spec == "spec":
        py_change_parameter.change_python_parameter(root, "Spectrum_cycles", "5")
        py_change_parameter.change_python_parameter(root, "Photons_per_cycle", "1e6")
    else:
        exit(1)

    start = time.time()
    command = "Setup_Py_Dir; mpirun -np {} {} {} {} >> output.out".format(ncores, py, flags, root)
    print(command)

    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = cmd.communicate()
    stdout = out.decode("utf-8")
    if verbose:
        print(stdout)
    stderr = out.decode("utf-8")
    if stderr:
        print(stderr)
    end = time.time()

    print("Total run time: {} s".format(end - start))

    if ion_or_spec == "ion":
        converged = check_convergence(root)
        ncycles = len(converged)
        print("Convergence of simulation:\n")
        for i in range(ncycles):
            print("{}/{}: converged fraction {}", i + 1, ncycles, converged[i])
        if converged == -1:
            print(ITS_A_MYSTERY)
            return 0
        elif converged[-1] >= CLIM:
            print(CONVERGED)
            return 0
        elif converged[-1] < CLIM:
            print(NOT_CONVERGED)
            for i in range(ncycles):
                if converged[i] >= CLIM:
                    break
            print("But.. simulation did converge on cycle {} with {}".format(i + 1, converged[i]))
            return 1

    return


def main():
    """
    Main control function.
    """

    p = argparse.ArgumentParser(description="Helper script to run simulations on Python on Iridis")
    p.add_argument("py", type=str, help="The path to the Python executable.")
    p.add_argument("flags", type=str, help="Any flags to pass to the Python executable.")
    p.add_argument("root", type=str, help="The root name of the Python simulation to run.")
    p.add_argument("ncores", type=str, help="The number of cores to run Python with.")
    p.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging.")
    args = p.parse_args()

    py = args.py
    flags = args.flags
    root = args.root
    ncores = args.ncores
    verbose = False  # args.verbose

    rc = run_python(py, flags, root, ncores, "ion", verbose)
    if rc == 0:  # If converged, do some spectral cycles
        run_python(py, flags, root, ncores, "spec", verbose)

    return


if __name__ == "__main__":
    main()
