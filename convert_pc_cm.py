#!/usr/bin/env python3

"""

Usage:
    convert_pc_cm.py 350 -M

    This will convert 350 Mpc into cm.

"""


import argparse
from sys import exit


verbose = False


def parse_distance():
    """
    Parse arguments from the command line

    Parameters
    ----------
    None

    Returns
    -------
    dist_pc: float
        The distance in parsecs to be converted.
    OoM_value: float
        The order of magnitude of the distance.
    """

    p = argparse.ArgumentParser(description="Convert a parsec distance to cm")
    p.add_argument("distance", type=float, help="The distance in parsecs")
    p.add_argument("-k", "-K", action="store_true", help="Use for kilo parsecs")
    p.add_argument("-m", "-M", action="store_true", help="Use for Mega parsecs")
    p.add_argument("-g", "-G", action="store_true", help="Use for Giga parsecs")
    p.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    args = p.parse_args()

    # Parse the distance in pc from the command line
    dist_pc = args.distance

    # Set up the order of magnitude factor
    OoM_value = 1  # default value
    if args.k:
        OoM_value = 1e3
    if args.m:
        OoM_value = 1e6
    if args.g:
        OoM_value = 1e9

    # Set up optional variables
    if args.verbose:
        global verbose
        verbose = True

    return dist_pc, OoM_value


def pc_to_cm(dist_pc, OoM):
    """
    Convert a distance in parsecs to cm.

    Parameters
    ----------
    dist_pc: float
        The distance in parsecs to be converted
    OoM: float
        The order of magnitude of the distance
    Returns
    -------
    dist_cm: float
        The distance in centimetres.
    """

    # Taken from https://www.unitconverters.net/length/parsec-to-centimeter.htm :^)
    pc_in_cm = 3.08567758128e18
    if verbose:
        print("{} pc in a cm".format(pc_in_cm))
        print("Value taken from https://www.unitconverters.net/length/parsec-to-centimeter.htm")
    dist_cm = dist_pc * OoM * pc_in_cm

    return dist_cm


def main():
    """
    Main steering function.
    """

    dist_pc, OoM = parse_distance()
    dist_cm = pc_to_cm(dist_pc, OoM)
    print("{} pc = {} cm".format(dist_pc * OoM, dist_cm))

    return


if __name__ == "__main__":
    main()
