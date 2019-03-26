#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Remove all data directories from Python simulations - this was written because Dropbox is able to follow symbolic
links on Linux >:(
"""


from sys import argv
from subprocess import Popen, PIPE


def remove_data_dir(search_dir:str="~/PySims")->int:
    """
    Search recursively from search_dir and downwards for directories, files and symlinks named data.

    Parameters
    ----------
    search_dir: str
        The base directory of which to search recursively for Python data symbolic links.

    Returns
    -------
    ndel: int
        The number of symbolic links which were deleted
    """

    # - type l will only search for symbolic links named data
    cmd = "cd {}; find . -type l -name 'data'".format(search_dir)
    stdout, stderr = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    stdout = stdout.decode("utf-8")
    stderr = stderr.decode("utf-8")

    if stderr:
        print("Message sent to stderr:")
        print(stderr)
    if stdout:
        print("Deleting data symbolic links in the following directories:\n\n{}".format(stdout))
    else:
        print("No data symlinks to delete")
        exit(0)

    # Create a hardcoded path and subprocess to remove each file
    ndel = 0
    dirs = stdout.split()
    for i in range(len(dirs)):
        dir = search_dir + dirs[i][1:]
        cmd = "rm {}".format(dir)
        stdout, stderr = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True).communicate()
        stdout = stdout.decode("utf-8")
        if stdout:
            print(stdout)
        stderr = stderr.decode("utf-8")
        if stderr:
            print(stderr)
        else:
            ndel += 1

    return ndel


if __name__ == "__main__":
    if len(argv) > 1:
        remove_data_dir(argv[1])
    else:
        remove_data_dir()
