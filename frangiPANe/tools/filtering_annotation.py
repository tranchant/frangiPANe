#!/usr/bin/env python3

# License GNU/GPL
# Intellectual property belongs to IRD
# written by Christine Tranchant-Dubreuil (UMR DIADE/itrop platform - IRD)

"""


    :Example:

    >>> python3 /home/christine/Script/Python/cdhitVsAnchoring.py


"""

import argparse
import sys
import os.path
import pprint as pp


from Bio import SeqIO

def main():

    parser = argparse.ArgumentParser(usage='%(prog)s [options]',
                                     description='Python script to generate fasta file ',
                                     epilog="""Author : C. Tranchant-Dubreuil (UMR DIADE / Rice Team - IRD) \n
                                                  ,--------------------""")

    parser.add_argument('-d', '--dir', help='Directory with gff and fasta merged', required=True)
    parser.add_argument('-i', '--id', help='List of contigs to select', required=True)
    args = parser.parse_args()

    maker_dir = args.dir
    id_file = args.id

    for file_name in os.listdir(maker_dir):
        print(file_name)
        if not "samtoolsFlagstat" in file_name:
            continue

        sample = os.path.basename(file_name).split('.')[0]
        newLine = f"\n{sample}"

        #with open(os.path.join(output_dir, file_name), "r") as stat:
        #    for line in stat:
        #        line = line.rstrip()


if __name__ == "__main__":
    main()