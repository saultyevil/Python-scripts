#!/usr/bin/env python3


"""
Solar mass fractions:
    X = 0.7381
    Y = 0.2485
    Z = 0.0134
    
We can just use,
    X = 0.74
    Y = 0.24
    Z = 0.02
for simplicities sake...
"""


import numpy as np
from sys import exit


verbose = True
extraVerbose = False


def readOpal(filename):
    """
    Read the Opal Opacity Tables into various lists and arrays.
    
    Parameters
    ----------
    filename: str
        The file path to the Opal Opacity table.
    
    Returns
    -------
    opalRMO: array of floats
        Blah
    opalMassFractions: array of floats
        Blah
    """

    print("Reading Opal Opacity Tables into memory")

    # Read in the Opal tables    
    try:    
        with open(filename, "r") as f:
            opal = f.readlines()
    except IOError:
        print("-----------------------------------------------------------------------")
        print("Couldn't find Opal Opacity table for file path {}".format(filename))
        print("-----------------------------------------------------------------------")
        exit(1)

    # Initialise various variables and data structures
    nLogT = 70      # Number of log Temperature points (nrows)
    nLogR = 19      # Number of log R points (ncols)
    nLines = 77     # Number of lines per table
    nTables = 126   # Number of tables
    tableId = 241   # Used to index from the first table

    # Data structures to hold Opal table
    opalRMO = np.zeros((nTables, nLogT+1, nLogR+1))
    opalMassFractions = np.zeros((nTables, 3))
    
    # Create array of the logR range
    logR = np.array(opal[tableId+4].split()[1:], dtype=float)
    
    nTable = 1
    for table in range(nTables):
        # Get the mass fractions
        for i in range(len(opal[tableId])):
            if opal[tableId][i] == "X" and opal[tableId][i+1] == "=":
                X = opal[tableId][i+2:i+8]
                if extraVerbose:    
                    print("X = {}".format(opal[tableId][i+2:i+8]))
            if opal[tableId][i] == "Y" and opal[tableId][i+1] == "=":
                Y = opal[tableId][i+2:i+8]
                if extraVerbose:
                   print("Y = {}".format(opal[tableId][i+2:i+8]))
            if opal[tableId][i] == "Z" and opal[tableId][i+1] == "=":
                Z = opal[tableId][i+2:i+8]
                if extraVerbose:
                    print("Z = {}".format(opal[tableId][i+2:i+8]))
        # X Y Z
        opalMassFractions[table, 0] = X
        opalMassFractions[table, 1] = Y
        opalMassFractions[table, 2] = Z
        
        # Assign the Opal data to numpy arrays
        opalRMO[table, 0, 1:] = logR
        for i in range(nLogT):
            line = opal[tableId+6+i].split()
            # Due to the Opal tables not being sqaure, keep appending 0 to the
            # end of the list to avoid dimension miss matches
            nAppends = 0
            while len(line) != nLogR+1:
                line.append(0)
                nAppends += 1
            if extraVerbose:
                print("Appended {} 0's to line {} in table {}".format(nAppends, i, nTable))
            opalRMO[table, 1+i, :] = line
        tableId += nLines
        nTable += 1

    print("Opal Opacity Tables loaded")

    return opalRMO, opalMassFractions


def readLowTempOpac(filenameOpac, filenameSets):
    """
    Read in the low temperature RMO data and set data into arrays.
    
    Parameters
    ----------
    filenameOpac: str
        The file path to the low temperature opacity data
    filenameSets: str
        The file path tot he low temperature opacity sets
    
    Returns
    -------
    lowTempRMO: array of floats
        Blah
    LowTempMassFractions: array of floats
        Blah
    """
    
    print("Reading LA08 Low Temperature Opacities into memory")
    
    # Read in the opacity data
    try:
        # LowTempHeaders = ["Zinit", "Set", "logT", "logR...."]
        with open(filenameOpac, "r"):
            lowTemp = np.loadtxt(filenameOpac)
    except IOError:
        print("-----------------------------------------------------------------------")
        print("Couldn't find Low Temperature Opacity table for file path {}".format(filenameOpac))
        print("-----------------------------------------------------------------------")
        exit(1)

    # Read in the "sets" data
    try:
        # LowTempSetsHeaders = ["Zinit", "Set", "X", "Y", "Z", "C/O", "12C", "14N", "alpha"]
        with open(filenameSets, "r"):
            lowTempSets = np.loadtxt(filenameSets)
    except IOError:
        print("-----------------------------------------------------------------------")
        print("Couldn't find Low Temperature Opacity sets table for file path{}".format(filenameSets))
        print("-----------------------------------------------------------------------")
        exit(1)
        
    # Generate the logT and logR spacing
    logR = np.arange(-7, 1.05, 0.5)
    
    # Reshape the array into the number of tables
    nLogT = 18
    lowTempRows = lowTemp.shape[0]
    lowTempCols = lowTemp.shape[1]
    nTables = int(lowTempRows / nLogT)
    lowTemp = np.reshape(lowTemp, (nTables, nLogT, lowTempCols))
    
    # Reshape the arrays into the format I want
    lowTempRMO = np.zeros((nTables, nLogT+1, lowTempCols-2))
    for table in range(nTables):
        lowTempRMO[table, 0, 1:] = logR
        lowTempRMO[table, 1:, :] = lowTemp[table, :, 2:]
    
    # Set X Y Z
    lowTempMassFractions = np.zeros((nTables, 4))
    for i in range(nTables):
        lowTempMassFractions[i, 0] = lowTempSets[i, 1]
        lowTempMassFractions[i, 1] = lowTempSets[i, 2]
        lowTempMassFractions[i, 2] = lowTempSets[i, 3]
        lowTempMassFractions[i, 3] = lowTempSets[i, 4]
    
    print("LA08 Low Temperature Opacities loaded")
    
    return lowTempRMO, lowTempMassFractions


def spliceTables(opal, opalMF, lowTemp, lowTempMF, X, Y, Z):
    return "Fuck"


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
    
    # Hardcode the file paths for now ;-)
    filenameOpal = "/home/saultyevil/Dropbox/Snake/examples/GN93hz"
    filenameOpac = "/home/saultyevil/Dropbox/Snake/LowTempOpac/opac.dat"
    filenameSets = "/home/saultyevil/Dropbox/Snake/LowTempOpac/set.dat"

    # Read in the tables into memory
    opalRMO, opalMassFrac = readOpal(filenameOpal)
    lowTempRMO, lowTempMassFrac = readLowTempOpac(filenameOpac, filenameSets)
    
    # Now splice the table together given the desired mass fractions
    X = 0.74
    Y = 0.24
    Z = 0.02
    newTable = spliceTables(opalRMO, opalMassFrac, lowTempRMO, lowTempMassFrac, X, Y, Z)

    return opalRMO, opalMassFrac, lowTempRMO, lowTempMassFrac, newTable


if __name__ == "__main__":
    opalRMO, opalMassFrac, lowTempRMO, lowTempMassFrac, newTable = main()
        