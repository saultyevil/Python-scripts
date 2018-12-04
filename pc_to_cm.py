#!/usr/bin/env python3

"""

Usage:
    pc_to_cm.py 350 -OoM M

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
    p.add_argument("-OoM", type=str, nargs="?", action="store", help="The prefix for the order of magnitude of the distance")
    p.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    args = p.parse_args()

    # Parse the distance in pc from the command line
    dist_pc = args.distance

    # Set up optional variables
    if args.verbose:
        global verbose
        verbose = True
    OoM_value = 1  # default value
    if args.OoM:
        OoM_str = args.OoM.lower()
        if OoM_str == "k":
            OoM_value = 1e3
        elif OoM_str == "m":
            OoM_value = 1e6
        elif OoM_str == "g":
            OoM_value = 1e9
        else:
            print("Only Kilo, Mega and Giga supported!")
            exit(1)

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

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    dist_pc, OoM = parse_distance()
    dist_cm = pc_to_cm(dist_pc, OoM)
    print("{} pc = {} cm".format(dist_pc * OoM, dist_cm))

    return


if __name__ == "__main__":
    main()
