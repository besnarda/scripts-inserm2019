#!/usr/bin/python3
# -*- coding: utf-8 -*-


#-------------------------------------
#
#  SCRIPT PYTHON frame calculator
#					
#  INPUT:
#	 - start of the gene
#  - start of your peak
#  OUTPUT:
#  -The frame you need to translate it correctly with transeq (EMBOSS)
#  					Author : Alban Besnard, albanbesnard@hotmail.fr
#-------------------------------------

try:
  import argparse
except ImportError:
  print("oops, the import /argparse/ didn't work")


parser = argparse.ArgumentParser()
parser.add_argument('-peak_start', required=True,type=int, help='position in bp where your peak start')
parser.add_argument('-gene_start', required=True,type=int, help='position in bp where your gene start')

args = parser.parse_args()
peak_start=args.peak_start
gene_start=args.gene_start

output = ((int(gene_start)-int(peak_start) ) % 3)+1
print(output)