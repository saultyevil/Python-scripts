#!/usr/bin/env python3

"""
Solar mass fractions:
    X = 0.7381
    Y = 0.2485
    Z = 0.0134
"""


import argparse
import numpy as np
from sys import exit
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt


# Global variables
verbose = False
show_comparison = False
splicingTemp = 3.8
opalTable = "opal_opac.dat"
la08Table = "la08_opac.dat"
la08TableSets = "la08_sets.dat"
newTableOutputName = "largerT_opacity.dat"
X = 0.7
Z = 0.02


def parse_inputs():
    """
    Parse various options from the command line
    """

    p = argparse.ArgumentParser(description="Splice two opacity tables together :-)")
    p.add_argument("-X", type=float, nargs="?", action="store", help="Hydrogen mass fraction")
    p.add_argument("-Z", type=float, nargs="?", action="store", help="Metals mass fraction")
    p.add_argument("-spliceT", type=float, nargs="?", action="store", help="The temperature to splice the tables between")
    p.add_argument("-opal_location", type=str, nargs="?", action="store", help="Location of the Opal opacity tables")
    p.add_argument("-low_temp_location", type=str, nargs="?", action="store", help="Location of the LA08 Low Temperature Opacity Tables")
    p.add_argument("-low_temp_sets", type=str, nargs="?", action="store", help="Location of the LA08 Low Temperature sets data")
    p.add_argument("-output_name", type=str, nargs="?", action="store", help="The name of the new opacity table")
    p.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    p.add_argument("-s", "--show_comparison", action="store_true", help="Show a plot of the log(RMO) against temperature")

    args = p.parse_args()

    if args.verbose:
        global verbose
        verbose = True
    if args.spliceT:
        global splicingTemp
        splicingTemp = args.spliceT
    if args.opal_location:
        global opalTable
        opalTable = args.opal_location
    if args.low_temp_location:
        global la08Table
        la08Table = args.low_temp_location
    if args.low_temp_sets:
        global la08TableSets
        la08TableSets = args.low_temp_sets
    if args.output_name:
        global newTableOutputName
        newTableOutputName = args.output_name
    if args.show_comparison:
        global show_comparison
        show_comparison = True
    if args.X:
        global X
        X = args.X
    if args.Z:
        global Z
        Z = args.Z

    return


def writeTable(table, outputName):
    """
    Write the new table out to disk.
    """

    print("- Writing table out to {}".format(outputName))

    with open(outputName, "w") as out:
        out.write("                                                                              logR\n")
        out.write("{:^7}\n".format("logT"))
        for i in range(table.shape[0]):
            for j in range(table.shape[1]):
                out.write("{:^+7.3f} ".format(table[i, j]))
            out.write("\n")

    return


def readOpal(filename):
    """
    Read in the Opal Opacity Tables
    """

    print("- Reading Opal Opacity")
    try:
        with open(filename, "r") as f:
            opal = f.readlines()
    except IOError:
        print("-----------------------------------------------------------------------")
        print("Couldn't find Opal Opacity table for file path {}".format(filename))
        print("-----------------------------------------------------------------------")
        exit(1)

    nLogT = 70      # Number of log Temperature points (nrows)
    nLogR = 19      # Number of log R points (ncols)
    nLines = 77     # Number of lines per table
    nTables = 126   # Number of tables
    tableId = 241   # Used to index from the first table

    # Arrays to hold Opal table
    opalRMO = np.zeros((nTables, nLogT+1, nLogR+1))
    opalMassFractions = np.zeros((nTables, 3))

    # Read in the table
    nTable = 1
    for table in range(nTables):
        # First, we'll get the mass fractions for each table
        X = Y = Z = 0
        for i in range(len(opal[tableId])):
            if opal[tableId][i] == "X" and opal[tableId][i+1] == "=":
                X = opal[tableId][i+2:i+8]
            if opal[tableId][i] == "Y" and opal[tableId][i+1] == "=":
                Y = opal[tableId][i+2:i+8]
            if opal[tableId][i] == "Z" and opal[tableId][i+1] == "=":
                Z = opal[tableId][i+2:i+8]

        opalMassFractions[table, 0] = X
        opalMassFractions[table, 1] = Y
        opalMassFractions[table, 2] = Z

        # Add the Opal data into the arrays
        opalRMO[table, 0, 1:] = np.array(opal[tableId+4].split()[1:], dtype=float)
        for i in range(nLogT):
            line = opal[tableId+6+i].split()
            # Due to the Opal tables not being square, keep appending 0 to the end of the list to avoid dimension
            # miss match
            nAppends = 0
            while len(line) != nLogR + 1:
                line.append(0)
                nAppends += 1
            if verbose and nAppends != 0:
                print("Appended {} 0's to line {} in table {}".format(nAppends, i, nTable))
            opalRMO[table, 1+i, :] = line
        tableId += nLines
        nTable += 1

    return opalRMO, opalMassFractions


def readLowTempOpac(filenameOpac, filenameSets):
    """
    Read in the low temperature RMO data
    """

    print("- Reading LA08 Low Temperature Opacity")
    try:
        with open(filenameOpac, "r"):
            lowTemp = np.loadtxt(filenameOpac)
    except IOError:
        print("-----------------------------------------------------------------------")
        print("Couldn't find Low Temperature Opacity table for file path {}".format(filenameOpac))
        print("-----------------------------------------------------------------------")
        exit(1)
    try:
        with open(filenameSets, "r"):
            lowTempSets = np.loadtxt(filenameSets)
    except IOError:
        print("-----------------------------------------------------------------------")
        print("Couldn't find Low Temperature Opacity sets table for file path{}".format(filenameSets))
        print("-----------------------------------------------------------------------")
        exit(1)

    nLogT = 18
    lowTempRows = lowTemp.shape[0]
    lowTempCols = lowTemp.shape[1]
    nTables = int(lowTempRows / nLogT)
    lowTemp = np.reshape(lowTemp, (nTables, nLogT, lowTempCols))
    lowTempMassFractions = lowTempSets[:, 1:5]

    # Reshape the arrays into the format I want
    logR = np.arange(-7, 1.05, 0.5)
    lowTempRMO = np.zeros((nTables, nLogT+1, lowTempCols-2))
    for table in range(nTables):
        lowTempRMO[table, 0, 1:] = logR
        lowTempRMO[table, 1:, :] = lowTemp[table, :, 2:]

    return lowTempRMO, lowTempMassFractions


def findLowTempIndex(lowTempMassFrac, X, Z):
    """
    Find the index for the LA08 data which has the same mass fractions X and Z.
    """

    idx = 0
    for i in range(lowTempMassFrac.shape[0]):
        if lowTempMassFrac[i, 1] == X and lowTempMassFrac[i, 3] == Z:
            idx = i
            break

    if idx == 0:
        print("\nUh oh! Couldn't find default LA08 table for X = {} and Z = {}".format(X, Z))
        print("4D interpolation hasn't been implemented because I'm being lazy")
        exit(1)

    return idx


def createNewTable(outputName, lowTempRMO, opalRMO, idx, X, Z):
    """
    Create the new opacity table featuring both low and high temperature RMO
    """

    print("- Generating new opacity table for X = {} and Z = {}\n".format(X, Z))

    #  Determine the total number of logT values in both tables
    nLowT = lowTempRMO.shape[1] - 1
    nOpalT = opalRMO.shape[1] - 1
    nLogT = nLowT + nOpalT

    # Read in the logT values form the table since Opal doesn't use a uniform spacing between logT
    logT = np.zeros(nLogT)
    logT[:nLowT] = lowTempRMO[0, 1:, 0]
    logT[nLowT:] = opalRMO[0, 1:, 0]

    # Find the unique values of logT and the index for the intercept between LA08 and Opal
    logT = np.unique(logT)
    nLogT = len(logT)
    splicingTemp = 3.8
    if 3.6 <= splicingTemp <= 3.9:
        x, = np.where(logT == splicingTemp)
        k = int(x)
    else:
        print("3.6 <= splicingTemp <= 3.9")
        exit (1)

    # Generate the logR values, there should hopefully be 17
    logR = np.arange(-7, 1.5, 0.5)
    nLogR = len(logR)

    # Create the new table and fill in the logR header and logT columns
    newTable = np.zeros((nLogT+1, nLogR+1))
    newTable[1:, 0] = logT
    newTable[0, 1:] = logR

    #
    # As I'm being lazy, we'll match using the default mass fractions in LA08 rather than writing 4D interpolation
    # routine to do this correctly. I may change this in the future to do it correctly...
    #

    # Add the LA08 data to the table.. not using vector operations as this *should* be temp code
    for i in range(k):
        for j in range(len(logR)):
            newTable[1+i, 1+j] = lowTempRMO[idx, 1+i, 1+j]

    #
    # I will generate the data for the LA08 mass fractions using the Opal interpolation function which was provided.
    # So.. for this you will need the Fortran code named opal_interp.f. It takes in 4 arguments, T6, R, X and Z, which
    # have the same meaning as in the Opal code. I call this program using subprocess for each grid point in the logT
    # logR table, so it will take a while as it has to read in the data table each time :^)
    #

    for i in range(nOpalT-1):
        T6 = 1e-6 * 10 ** logT[i+k]
        for j in range(nLogR):
            R = 10 ** logR[j]
            # TODO: error checking for opal to make sure it can find the opacity tables
            callOpal = "libs/opal {:e} {:e} {} {}".format(T6, R, X, Z)
            stdout, stderr = Popen(callOpal, stdout=PIPE, stderr=PIPE, shell=True).communicate()
            try:
                opalRMO = float(stdout.decode("utf-8"))
            except ValueError:
                if verbose:
                    print("logT = {:1.2e} or logR = {:1.2e} out of table range".format(logT[i+k], logR[j]))
                    print("Setting to 9.999 for table element [{}, {}]".format(i+1+k, 1+j))
                opalRMO = 9.999
            newTable[1+i+k, 1+j] = opalRMO
        if i % 5 == 0 and i != 0:
            print("\t- Row {} of {} completed".format(i, nLogT-k))

    print("\t- Table completed\n")
    writeTable(newTable, outputName)

    return newTable


def plot_comparisons(newTable, la08Table):
    """
    Plot comparisons of the new table
    """

    print("- Creating comparison plots as requested")

    nrows = 2
    ncols = 2
    cols = [1, 7, 12, 17]
    titles = ["logR = -7.0", "logR = -4.0", "logR = -1.5", "logR = 1.0"]
    fig, ax = plt.subplots(nrows, ncols, figsize=(11.69, 8.27))

    k = 0
    for i in range(nrows):
        for j in range(ncols):
            colIdx = cols[k]
            title = titles[k]
            ax[i, j].plot(newTable[:13, 0], newTable[:13, colIdx], "-x", label="New table")
            ax[i, j].plot(la08Table[2:15, 0], la08Table[2:15, colIdx], "--+", label="LA08")
            ax[i, j].set_title(title)
            ax[i, j].set_xlabel("log T [K]")
            ax[i, j].set_ylabel("log $\kappa_{R}$ [$cm^{2}g^{-1}$]")
            ax[i, j].set_xlim(newTable[0, 0], newTable[12, 0])
            ax[i, j].legend()
            k += 1

    fig.tight_layout()
    plt.show()

    return


def main():
    """
    Main function
    """

    print("Welcome to the degenerate code Jack 'Oliver' Cuthbertson\n")
    print("We will be splicing together two opacity tables together to create")
    print("one large temperature range table\n")

    parse_inputs()

    # Read in the tables into file. I assume that the opacity tables are within the current directory or within a
    # sub-directory:
    #  - The directory paths can be changed with some switches on the command line
    opalRMO, opalMassFrac = readOpal(opalTable)
    lowTempRMO, lowTempMassFrac = readLowTempOpac(la08Table, la08TableSets)

    # Now splice the table together given the desired mass fractions
    idx = findLowTempIndex(lowTempMassFrac, X, Z)
    newTable = createNewTable(newTableOutputName, lowTempRMO, opalRMO, idx, X, Z)
    if show_comparison:
        plot_comparisons(newTable[3:, :], lowTempRMO[idx, 1:, :])

    print("\nEnjoy!")

    return

if __name__ == "__main__":
    main()
