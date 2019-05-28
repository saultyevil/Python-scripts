#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import py_plot_util
import py_change_parameter


def update_parameter_file(pfs, parameter_name, parameter_value):
    """

    Parameters
    ----------
    pfs
    parameter_name
    parameter_value

    Returns
    -------

    """


    if type(pfs) != list:
        print("update_parameter_file: type of pfs is not a list")
        return

    assert(type(parameter_name) is str)
    assert(type(parameter_value) is str)

    for pf in pfs:
        py_change_parameter.change_python_parameter(pf, parameter_name, parameter_value)

    return


if __name__ == "__main__":
    pfs = list(py_plot_util.find_pf())
    parameter = "Spectrum_cycles"
    value = "0"
    update_parameter_file(pfs, parameter, value)
