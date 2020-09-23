#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:16:21 2020

Python script to convert SAM file to BAM

Problems: 
    -POS is 1-based 
    -not mapped on reference genome
    
What we want: 
    -chr - 3rd column 
    -start pos (col 4, but 0 based)
    -end pos (col 3 _ len(col10))
    -name (col 1)
    -score (col 5)
    -strand (.) (Can put . if don't know format')

@author: Sophie
"""

#path for input samfile and bedfile output
#samfilepath = '/Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/ERR1755082.test.sam'
#bedfilepath = '/Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/ERR1755082.test.bed'

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', dest='samfilepath', help='input file path')
parser.add_argument('--output', '-o', dest='bedfilepath', help='output file path')
args=parser.parse_args()



with open(args.samfilepath, 'r') as samfile: # open the input file (defined above)
    with open(args.bedfilepath, 'w') as bedfile: # create the output file
        for line in samfile: # read the input file line by line
            if line[0] == '@': # headers in SAM starts with @, so we skip it
                pass
            else:
                col = line.split() # parse the line into a list of fields / columns
                #     - chr (column 3)
                chrom = col[2]
                #     - start pos (column 4 but need to convert from 1-based SAM file to 0-based BED file)
                startpos = int(col[3]) - 1
                #     - end pos (column 4 + len(column 10))
                endpos = int(col[3]) + len(col[9])
                #     - name (column 1)
                name = col[0]
                #     - score (column 5)
                score = col[4]
                #     - strand (.)
                strand = '.'
                bedfile.write(f'{chrom}\t{startpos}\t{endpos}\t{name}\t{score}\t{strand}\n')
        
   
   #on command line
   #(obds-py3) Sophies-MacBook-Pro-4:Python Sophie$ python samtobed.py -i ../../ERR1755082.test.sam -o ../../ERR1755082.test.bed  
   
   