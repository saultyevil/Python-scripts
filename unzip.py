#!/usr/bin/env python3

import argparse
from subprocess import Popen, PIPE

"""
Unzip a bunch of .zip files.. 

Usage:
    ./unzip.py file [--verbose]
    python unzip.py file [--verbose]
    
I suggest using ls > file to get all of the zips in a directory
"""


def unzip(filename, verbose):
    """
    Uses subprocess to communicate with the shell. Uses the Linux unzip command to unzip zip files into a directory
    of the same name

    Parameters
    ----------
    filename: str
        The file path of the file containing the zip files to unzip.
    verbose: bool
        If True, extra output will be printed

    Returns
    -------
    None
    """

    try:
        f = open(filename, "r")
        files = f.readlines()
        f.close()
    except IOError:
        print("Couldn't open file {}".format(filename))
        return -1
    nfiles = len(files)
    print("\n--------------------------------------------------------------------------------")
    print("Unzipping {} files:\n".format(nfiles))
    for i in range(nfiles):
        print("\t - {}".format(files[i].strip()))
    print("\n--------------------------------------------------------------------------------\n")
    for i in range(nfiles):
        file = files[i].strip()
        for j in range(len(file)):
            if file[j] == ".":
                break
        dirname = file[0:j]
        options = "-d {}".format(dirname)
        unzip = "unzip {} {}".format(file, options)
        print("{}/{}: unzipping {}".format(i + 1, nfiles, file, dirname))
        stdout, stderr = Popen(unzip, stdout=PIPE, stderr=PIPE, shell=True).communicate()
        stdout = stdout.decode("utf-8")
        if verbose:
            print(stdout)
    print("Done!")
    print("--------------------------------------------------------------------------------")
    return

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str, 
                        help="A file containing containing the pathes for files to unzip")
    parser.add_argument("--verbose", help="Increase verbosity, i.e. output stdout",
                        action="store_true")
    args = parser.parse_args()
    name = args.filename
    unzip(name, args.verbose)
