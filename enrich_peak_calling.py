#!/usr/bin/python3
# -*- coding: utf-8 -*-


#-------------------------------------
#
#  SCRIPT PYTHON enrich_peak_calling
#
#  From a gff file and the output of MACS 2.1 (peak calling)
#  Create an enriched file with the name of gene
#  Possible addon: use a file with more info on genes to add on the enriched file!
#
#  INPUT:
#	 - .gff file contaning the reference
#    - .xls file output from MACS
#    - OPTIONAL a supplemenntary .txt file with info on genes
#
#  					Author : Alban Besnard, albanbesnard@hotmail.fr
#-------------------------------------



#------------------------------------------------------------------------------------------------------
### STEP 0 : Importation, Description and Verification of input/output
#------------------------------------------------------------------------------------------------------

import os
import sys
import re
import gffutils


try:
  import argparse
except ImportError:
  print("oops, the import /argparse/ didn't work")


My_description="\nFrom a .xls file of MACS 2.1 peak calling and the .gff file of the reference, \nthis script will determine if each peak is within a gene. \nWhen it is available, it will add information about this gene \n output is a .xls MACS 2.1 like with some additional columns"

parser = argparse.ArgumentParser(description= My_description)
parser.add_argument('-ref', required=True, help='Reference genome with annotation used to make the peak calling (.gff format expected)')
parser.add_argument('-peak', required=True, help='File containing peak calling results (tabulation format expected)')
parser.add_argument('-supp', required=False, help='additional information about genes (tabulation format expected)')
parser.add_argument('-out', required=True, help='output file')


args = parser.parse_args()
fic_ref=args.ref
fic_peak=args.peak
fic_supp=args.supp
fic_out=args.out


#------------------------------------------------------------------------------------------------------
### STEP 1 : First tests
#------------------------------------------------------------------------------------------------------

# opening our peak file
peak = open(fic_peak,'r');
out = open(fic_out,'w');

# creating database for gffutils
#db = gffutils.create_db(fic_ref, dbfn='test.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
db = gffutils.FeatureDB('test.db', keep_order=True)

# we make a loop on each line of the peak file.
for line in peak:
    # All line begining with # are skipped and directly written in the output file.
    if line.startswith('#'):
        out.write(line)
        continue
    # there we need to say what columns we add (for the moment only the gene name)
    elif line.startswith('chr'):
        out.write(line.rstrip('\n') + '\tgene\n' )
        continue
    else
        # Here we enter the Core of our script
        linelist = re.split(r'\t+', line)
        chro = linelist[0]+".1"
        if chro == "Chromosome.1":
            chro = "CP000325.1"
        start = int(linelist[1])
        end = int(linelist[2])
        mid = int(linelist[4])
        print(mid)
        # find if the mid point is inside a gene
        for feature in list(db.features_of_type("gene")):
            if mid > feature.start and mid < feature.end and chro == feature.seqid:
                print(feature.id + feature.strand)

                continue



peak.close();
out.close();
