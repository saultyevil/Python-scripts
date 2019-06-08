#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import py_plot_util
import py_change_parameter


def update_parameter_file(pfs, parameter_name, parameter_value):
    """
    Uses py_change_parameter to update a bunch of provided parameter files.

    Parameters
    ----------
    pfs                    list[str]
                           A list containing the directories to the parameter
                           files to update
    parameter_name         str
                           The name of the parameter file to change
    parameter_value        str
                           The new value of the parameter
    Returns
    -------
    None

    """

    if type(pfs) != list:
        print("update_parameter_file: type of pfs is not a list")
        return

    assert(type(parameter_name) is str)
    assert(type(parameter_value) is str)

    for pf in pfs:
        print("Changing {} to {} in {}".format(parameter_name, parameter_value, pf))
        py_change_parameter.change_python_parameter(pf, parameter_name, parameter_value)

    return


if __name__ == "__main__":

    pfs = list(py_plot_util.find_pf_files())

    parameter = "Spectrum_cycles"
    value = "0"

    print("ENSURE THAT THE SCRIPT HAS BEEN EDITED APPROPRIATELY BEFORE RUNNING")
    input("Press a enter to continue...")

    update_parameter_file(pfs, parameter, value)
