#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
from PyPython import WindPlot
from PyPython import PythonUtils as Utils
from subprocess import Popen, PIPE


def plot_wind_saves():
    """
    Main function for the purpose that maybe this will one day want to be used
    in another script.
    """

    vars = ["t_e", "t_r", "ne", "ntot", "c4", "ip", "converge", "converging"]
    var_types = ["wind"] * len(vars)
    projection = "rectilinear"

    files = sorted(glob.glob("./python*.wind_save"), key=str.lower)
    nfiles = len(files)
    for i in range(nfiles):
        root = files[i][:files[i].find(".wind_save")]
        if root[:2] == "./":
            root = root[2:]
        print(root)
        sh = "Setup_Py_Dir; windsave2table -d {}".format(root)
        stdout, stderr = Popen(sh, stdout=PIPE, stderr=PIPE, shell=True).communicate()
        if stderr:
            print(stderr.decode("utf-8"))
            break
        # print(stdout.decode("utf-8"))
        input_file = "{}.0.master.txt".format(root)
        # WindPlot.rectilinear_wind(root, root, vars, var_types, "./", projection=projection, input_file=input_file,
        #                   verbose=True, plot_indices=True)

    Utils.remove_data_sym_links("./", verbose=True)

    return


if __name__ == "__main__":
    plot_wind_saves()
