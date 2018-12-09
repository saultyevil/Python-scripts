#!/usr/bin/env python3

"""
Solar mass fractions:
    X = 0.7381
    Y = 0.2485
    Z = 0.0134
"""


import numpy as np
from sys import exit
from subprocess import Popen, PIPE


verbose = False


def parse_inputs():
    return


def find_opacity_tables():
    """
    Find the file paths to the data which we need, hard coded for now
    """

    opalTable = "GN93hz"
    la08Table = "opac.dat"
    la08TableSets = "set.dat"

    opalFile = "/home/saultyevil/Dropbox/Snake/examples/GN93hz"
    la08File = "/home/saultyevil/Dropbox/Snake/LowTempOpac/opac.dat"
    la08Sets = "/home/saultyevil/Dropbox/Snake/LowTempOpac/set.dat"

    return opalFile, la08File, la08Sets


def writeTable(newTable, outputName):
    """
    Write the new table out to disk.
    """

    print("\t-Writing table out to {}".format(outputName))

    with open(outputName, "w") as out:
        out.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t  logR\n")
        out.write("{:^7}\n".format("logT"))
        for i in range(newTable.shape[0]):
            for j in range(newTable.shape[1]):
                out.write("{:^7.3f} ".format(newTable[i, j]))
            out.write("\n")

    return


def readOpal(filename):
    """
    Read in the Opal Opacity Tables
    """

    print("\t-Reading Opal Opacity")
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
            if verbose:
                print("Appended {} 0's to line {} in table {}".format(nAppends, i, nTable))
            opalRMO[table, 1+i, :] = line
        tableId += nLines
        nTable += 1

    return opalRMO, opalMassFractions


def readLowTempOpac(filenameOpac, filenameSets):
    """
    Read in the low temperature RMO data
    """

    print("\t-Reading LA08 Low Temperature Opacity")
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


def createNewTable(outputName, lowTempRMO, opalRMO, idx, X, Z):
    """
    Create the new opacity table featuring both low and high temperature RMO
    """

    print("\t-Generating new opacity table")

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
        k, = int(np.where(logT == splicingTemp))
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
        T6 = 10 ** logT[i+k] * 1e-6
        for j in range(nLogR):
            R = 10 ** logR[j]
            callOpal = "./opal {:e} {:e} {} {}".format(T6, R, X, Z)
            stdout, stderr = Popen(callOpal, stdout=PIPE, stderr=PIPE, shell=True).communicate()
            try:
                opalRMO = float(stdout.decode("utf-8"))
            except ValueError:
                if verbose:
                    print("logT = {:1.2e} or logR = {:1.2e} out of table range".format(lT, lR))
                    print("Setting to 9.999 for table element [{}, {}]".format(i+1+k, 1+j))
                opalRMO = 9.999
            newTable[1+i+k, 1+j] = opalRMO
        if i % 5 == 0:
            print("\t    - Row {} of {} completed".format(i+1, nLogT-k))

    print("\t    - Table completed")
    writeTable(newTable, outputName)

    return


def main():
    """
    Main function
    """

    print("Welcome to the degenerate code Jack Oliver Cuthbertson:")
    print("We will be splicing together two opacity tables for a large temperature range")

    # Read in the tables into file. I assume that the opacity tables are within the current directory or within a
    # sub-directory
    opalFile, la08File, la08Sets = find_opacity_tables()
    opalRMO, opalMassFrac = readOpal(opalFile)
    lowTempRMO, lowTempMassFrac = readLowTempOpac(la08File, la08Sets)

    # Now splice the table together given the desired mass fractions, which is hardcoded for now
    X = 0.70
    Z = 0.02
    idx = 70
    outputName = "final_table.txt"
    createNewTable(outputName, lowTempRMO, opalRMO, idx, X, Z)

    return

if __name__ == "__main__":
    main()
