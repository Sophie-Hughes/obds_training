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

import gzip 
import argparse
parser = argparse.ArgumentParser(description="script to convert sam to bed format")
parser.add_argument('--input', '-i', dest='samfilepath', help='input file path')
parser.add_argument('--output', '-o', dest='bedfilepath', help='output file path')
 # Provide a command line argument to pad the intervals in the bed file  
parser.add_argument('--padding', '-p', dest='bedfilepad', type=int, default=0, help='number of bases added to start and end of interval')
parser.add_argument('--fragment', '-f', dest='fragment', action="store_true", default=False, help='write out fragments rather than reads')
args=parser.parse_args() 





with open(args.samfilepath, 'r') as samfile: # open the input file (defined above)
    #Write out a gzip compressed file (.bed.gz) - 'wt' as it's a txt file
    with gzip.open(args.bedfilepath, 'wt') as bedfile: # create the output file
        for line in samfile: # read the input file line by line
            if line[0] == '@': # headers in SAM starts with @, so we skip it
                pass
            else:
                col = line.split() # parse the line into a list of fields / columns
                #     - chr (column 3)
                chrom = col[2]
                # if chrom doesn't have a '*', continue on. can also write if chrom == '*':
                    #continue --> this is also correct 
                if chrom != '*':
                         #     - score (column 5)
                   score = int(col[4])
                   if score > 0:
                    #     - start pos (column 4 but need to convert from 1-based SAM file to 0-based BED file)
                        if args.bedfilepad == 0:    
                            startpos = int(col[3]) - 1
                            endpos = int(col[3])+ len(col[9])
                        elif args.bedfilepad > 0:
                            startpos = int(col[3]) - (1 + args.bedfilepad) 
                         #     - end pos (column 4 + len(column 10))
                            endpos = int(col[3]) + len(col[9]) + args.bedfilepad
                #     - name (column 1)
                        name = col[0]
                #     - strand (.)
                        strand = '.'
                        #check fragment mode
                        if args.fragment:   #don't have to put = True because boolean will already evaluate this 
                            #tlen to avoid the re-pairs twice 
                            tlen = int(col[8])
                            if tlen > 0: 
                                endpos = int(col[3]) + tlen #if tlen greater than 0 we change the endpos
                                bedfile.write(f'{chrom}\t{startpos}\t{endpos}\t{name}\t{score}\t{strand}\n')
                        else: 
                            bedfile.write(f'{chrom}\t{startpos}\t{endpos}\t{name}\t{score}\t{strand}\n')
                                
                                
                 
                    
      
    
  
   
    
   
  
   