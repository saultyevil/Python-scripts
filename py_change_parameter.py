#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is a quick and simple script to change the value of a parameter in a Python
.pf file.

The script expects 3 arguments, as documented below.

Usage
   $ python py_change_parameter.py pf parameter value

        - pf            The name of the parameter file to edit
        - parameter     The name of the parameter to edit
        - value         The new value for the parameter
"""


from sys import argv, exit
from shutil import copyfile


def change_python_parameter(pf: str, parameter: str, value: str, verbose: bool = False, bakup: bool = True):
    """
    Search a parameter file for a given parameter and replaces the current value
    with a new value. This script will change the parameter file, even if the
    old and new parameter values are the same :-).

    Parameters
    ----------
    pf              str
                    The name of the parameter file to edit
    parameter       str
                    The name of the parameter to be edited
    value           str
                    The new value of the parameter
    verbose         bool, optional
                    Enable verbose logging
    bakup           bool, optional
                    If True, save a back up of the parameter file prior to edit

    Returns
    -------
    Returns a non-zero integer on non-successful exit.
    """

    assert(type(pf) == str)
    assert(type(parameter) == str)
    assert(type(value) == str)

    if pf.find(".pf") == -1:
        pf += ".pf"

    old = ""
    new = ""

    # Create back up file, in case things go to shit
    if bakup:
        copyfile(pf, pf + ".bak")

    # Open file, search for parameter and replace
    with open(pf, "r") as f:
        lines = f.readlines()
    nlines = len(lines)
    for l in range(nlines):
        if lines[l].find(parameter) != -1:
            old = lines[l]
            new = "{} {}\n".format(parameter, value)
            lines[l] = new
            break
    if old and new:
        if verbose:
            print("Changed parameter {} to value {}".format(parameter, value))
            print("OLD: {}".format(old.replace("\n", "")))
            print("NEW: {}".format(new.replace("\n", "")))
    else:
        print("Could not find parameter {} in {}".format(parameter, pf))
        exit(1)

    # Now write out modified lines to file
    with open(pf, "w") as f:
        f.writelines(lines)

    return


if __name__ == "__main__":
    if len(argv) != 4:
        print(__doc__)
        exit(0)
    else:
        root = argv[1]
        parameter = argv[2]
        value = argv[3]
    change_python_parameter(root, parameter, value, True)
