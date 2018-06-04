"""
This script transfers basic fusion granules between two globus endpoints.
Given a source directory and orbit, the script infers the location
of the BF granule based on a set of rules:

1. That the source directory contains subdirectories of the format:
   yyyy.mm
2. That each BF granule is named according to the BF naming convention
   and that it exists.

The script also assumes that the orbits you request are available in
the provided source directory.
"""

import argparse
import globuslite
import bfutils
import os

def main():
    parser = argparse.ArgumentParser(description='This script provides \
    an easy way to transfer Basic Fusion granules between two Globus \
    endpoints.', 
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('src_id', help='Source Globus endpoint ID',
        type=str)
    parser.add_argument('dest_id', help='Destination Globus endpoint ID',
        type=str )
    parser.add_argument('src_dir', help='Source path to the Basic Fusion \
        files. This directory must contain the the BasicFusion \
        subdirectories of the convention: "yyyy.mm" where yyyy is the \
        year, and mm is the month.', type=str )
    parser.add_argument('dest_dir', help="Destination directory. \
    Globus will mkdir this directory if it doesn't exist. Will contain \
    the yyyy.mm directories.", type=str)
    parser.add_argument('orbit_lower', help='Specify lower bound of orbit \
        range to transfer (inclusive).', type=int )
    parser.add_argument('orbit_upper', help='Specify upper bound of orbit \
        range to transfer (inclusive).', type=int )

    args = parser.parse_args()

    transfer = globuslite.Transfer( args.src_id, args.dest_id )

    for orbit in range( args.orbit_lower, args.orbit_upper + 1 ):
        stime = bfutils.file.orbit_start( orbit )

        bf_dir = stime[0:4] + '.' + stime[4:6]
        fname = bfutils.file.get_bf_filename(orbit)

        src_path = os.path.join( args.src_dir, bf_dir, fname )
        dest_path = os.path.join( args.dest_dir, bf_dir, fname )

        transfer.add_item( src_path, dest_path )

    transfer.submit()
        
if __name__ == '__main__':
    main()
