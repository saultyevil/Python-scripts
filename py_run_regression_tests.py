#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shutil
from subprocess import Popen, PIPE
import py_run_util


def get_latest_python():
    """
    TODO: need a better way to check for an update on the remote repo :-)
    """

    try:
        py_dir = os.environ["PYTHON_DEV"]
    except KeyError:
        print("PYTHON_DEV environment variable not set")
        sys.exit(1)

    print("Python test directory {}".format(py_dir))
    print("Pulling latest version of Python and compiling\n")

    pull = "cd {}/source; git reset --hard; git checkout dev; git pull; make all; make clean".format(py_dir)
    cmd = Popen(pull, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    if cmd.returncode:
        print("An error occured when trying to update and compile Python")
        print("error code {}".format(cmd.returncode))
        print(stderr.decode("utf-8"))
        sys.exit(1)

    return py_dir


def main():
    """

    """

    print("Running Python regression test")
    print("------------------------------\n")

    rc = shutil.which("regression.py")
    if not rc:
        print("regression.py not in $PATH")
        sys.exit(1)
    rc = shutil.which("regression_check.py")
    if not rc:
        print("regression_check.py not in $PATH")
        sys.exit(1)
    rc = shutil.which("regression_plot.py")
    if not rc:
        print("regression_plot.py not in $PATH")
        sys.exit(1)

    # Determine the average load in the last 15 minutes using uptime
    uptime = "uptime"
    cmd = Popen(uptime, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    stdout = stdout.decode("utf-8").split()
    try:
        load_avg = float(stdout[-1])
    except ValueError:
        print("Could not determine load average for last 15 minutes")
        sys.exit(1)

    if load_avg > 10:
        print("Load average > 10%, therefore the computer is probably doing something")
        sys.exit(1)

    try:
        regress_dir = os.environ["REGRESS_DIR"]
        print("Regression directory {}\n".format(regress_dir))
    except KeyError:
        print("REGRESS_DIR environment variable not set")
        sys.exit(1)

    # Pull in the latest version of Python and compile
    py_dir = get_latest_python()

    # Now we can run the regression script
    print("Running regression tests...\n")

    def_cores = 3
    py_ver = py_dir + "/bin/py"
    n_cores = py_run_util.get_num_procs(def_cores)
    regress = "cd {}; regression.py {} -np {}\n".format(regress_dir, py_ver, n_cores[1])
    print(regress)
    cmd = Popen(regress, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    stdout = stdout.decode("utf-8")
    stderr = stderr.decode("utf-8")

    print("Regression test return code: ", cmd.returncode)

    print(stdout)
    if stderr:
        print(stderr)

    return


if __name__ == "__main__":
    main()