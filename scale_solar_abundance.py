#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Used to scale some non-solar abundances to the weird scaling that astronomers
use where abundances are scaled wrt solar hydrogen in the sun or whatever.
"""

import numpy as np

# CNO Processed star
#   He ~ 2xsolar
#   C ~ 0.5xsolar
#   N ~ 7xsolar
#   O ~ 1xsolar


def main():
    """
    Given some input, scale...
    """

    origA = [10.99, 8.56, 8.05, 8.93]
    scale = [2.0, 0.5, 7.0, 1.0]
    scaledA = np.zeros(len(scale))

    assert(len(origA) == len(scale))

    for i in range(len(origA)):
        Nel = np.exp(origA[i] - 12) * 1e12
        Nel *= scale[i]
        scaledA[i] = np.log(Nel / 1e12) + 12

    return scaledA


if __name__ == "__main__":
    scaledA = main()
    print("The scaled abundances are: ", scaledA)