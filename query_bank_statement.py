#!/usr/bin/env python3

"""
The purpose of this script is to simply try and track what my housemates have
paid me for rent and utilities. Reads in a bank statement from CSV and searches
for transactions with the following names:
    - GREENSLADE
    - GUTHRIE
    - DEGREGORI
The script will then make a list of the transactions to try and simplify 
tracking what has and has not been paid.
Note that this script is expecting the statement to be a CSV file.
"""

from sys import exit
import numpy as np
import pandas as pd
import argparse as ap

COL_LEN = 126
VERBOSE = False
NAME_HEADER = "Transaction Description"
CREDIT_HEADER = "Credit Amount"
DEBIT_HEADER = "Debit Amount"

TOM_NAME = "GREENSLADE"
BEN_NAME = "GUTHRIE"
ENRICO_NAME = "DEGREGORI"
NAMES = [TOM_NAME, BEN_NAME, ENRICO_NAME]


def get_filename() -> str:
    """Parse the filename of the bank statement from the command line"""

    p = ap.ArgumentParser(description="Simple script to identify transactions made by three people in the provided"
                          " statement.")
    p.add_argument("filename", type=str, action="store", help="The filename of the bank statement.")
    p.add_argument("--name", type=str, action="store", help="Specify the name of a single person/company")
    p.add_argument("--nameheader", type=str, action="store", help="Define the name of the name column in the"
                   " statement.")
    p.add_argument("--creditheader", type=str, action="store", help="Define the name of the credit column in the" 
                   " statement.")
    p.add_argument("--debitheader", type=str, action="store", help="Define the name of the debit column in the"
                   "statement.")
    p.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging.")
    args = p.parse_args()

    if args.verbose:
        global VERBOSE
        VERBOSE = True
    if args.name:
        global NAMES
        NAMES = [args.name]
    if args.nameheader:
        global NAME_HEADER
        NAME_HEADER = args.nameheader
    if args.creditheader:
        global CREDIT_HEADER
        CREDIT_HEADER = args.creditheader
    if args.debitheader:
        global DEBIT_HEADER
        DEBIT_HEADER = args.debitheader

    return args.filename


def read_in_statement(filename: str) -> pd.DataFrame:
    """Read in a bank statement into a Pandas DataFrame"""

    try:
        df = pd.read_csv(filename)
    except IOError:
        print("Could not open the file {}. Does it exist?".format(filename))
        exit(1)

    return df


def get_transactions(df: pd.DataFrame, name: str) -> pd.DataFrame:
    """Get the transactions from a DataFrame for a specific person"""

    if name not in NAMES:
        print("{} is not an expected name, good luck!".format(name))

    try:
        sdf = df[df[NAME_HEADER].str.contains(name)]
    except KeyError:
        print("Unable to find entries for {} in column named {}. Do both exist?".format(name, NAME_HEADER))
        return pd.DataFrame()

    return sdf


def main():
    """The main function of the script"""

    pd.set_option('display.expand_frame_repr', False)

    filename = get_filename()
    df = read_in_statement(filename)
    
    print("-" * COL_LEN)
    print("Statement:", filename)
    print("Total transactions:", df.shape[0])

    hmates = []
    for name in NAMES:
        hmates.append(get_transactions(df, name))
    tot_credit = 0
    tot_debit = 0
    for i, trans in enumerate(hmates):
        credit_sum = np.sum(trans[CREDIT_HEADER])
        debit_sum = np.sum(trans[DEBIT_HEADER])
        tot_credit += credit_sum 
        tot_debit += debit_sum
        print("-" * COL_LEN)
        print("Transactions from:", NAMES[i])
        print(trans)
        print("Total credit: £{:4.2f}".format(credit_sum))
        print("Total debit: £{:4.2f}".format(debit_sum))
    print("-" * COL_LEN)
    print("Total credit: £{:4.2f}".format(tot_credit))
    print("Total debit: £{:4.2f}".format(tot_debit))
    print("-" * COL_LEN)

    return


if __name__ == "__main__":
    main()
