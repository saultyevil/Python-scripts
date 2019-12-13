#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xraydb
from consts import HEV
from typing import List

C = 299792458
ANGSTROM = 1e-10
COL_LEN = 80


def get_xray_edges(elements: List[str], wmin: float, wmax: float):
    """
    Using xraydb, return the absorbtion edges

    Parameters
    ----------
    elements: List[str]
        A list of the element symbols from which to query absorption edges.
    wmin: float
        The smallest wavelength edge to return
    wmax: float
        The largest wavelength edge to return

    Returns
    -------
    output_table: List[str]
        A table containing absorption edges.
        - Elem: the element
        - Energy: the photoionisation energy
        - Frequency: the frequency of the absorption edge
        - Wavelength: the wavelength of the absorption edge
    """

    element_absortion_edges_dicts = []
    for element in elements:
        edges = xraydb.xray_edges(element)
        element_absortion_edges_dicts.append(edges)

    output_table = []
    output_table.append("Elem {:15s} {:15s} {:15s}\n".format("Energy eV", "Frequency Hz", "Wavelength AA"))

    for i, edges in enumerate(element_absortion_edges_dicts):
        print("-" * COL_LEN)
        print("{}: \n".format(elements[i]))
        print("{:15s} {:15s} {:15s}".format("Energy eV", "Frequency Hz", "Wavelength AA"))
        keys = edges.keys()
        prev_key = "K"
        for key in keys:
            # This bit will skip edges which have the same energy, I hope
            if prev_key != key:
                if edges[prev_key][0] == edges[key][0]:
                    continue
            prev_key = key
            energy = edges[key][0]
            frequency = energy / HEV
            wavelength = C / frequency / ANGSTROM
            print("{:9.1f} {:1.12e} {:13.1f}".format(energy, frequency, wavelength))
            if wmin < wavelength < wmax:
                output_table_line = "{:4s} {:9.1f} {:1.12e} {:13.1f}\n".format(
                        elements[i], energy, frequency, wavelength)
                output_table.append(output_table_line)
        print()
    print("-" * COL_LEN)

    with open("xray_edges.txt", "w") as f:
        f.writelines(output_table)

    return output_table


# Wavelength limits in Angstroms
wmin = 100
wmax = 1000

# These are the current elements which we use in Python
elements = ["H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al", "Si", "S", "Ar", "Ca", "Fe"]
edges = get_xray_edges(elements, wmin, wmax)
