#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:29:00 2020

1)import pysam 
2)
3)read in SAM file in pysam 
4)iterate over SAM file line by line 
5)create bed file 
6) mandaroty fields for bed file: 
    1)chromosome
    2)cheom start (0) based
    3)Chrom End
    4)name
    5)score (mapping_quality)
    6)strand
    
    
4)name what BAM = 
3)remove header with pysam 
4)then iterate over every line 
5)
6)
7)

@author: Sophie
"""
#import pysam 
import pysam 
import argparse
#input filepath 
# bamfilepath = '/Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/obds_training/Python/ERR1755082.test.sort.bam'
# '/ifs/obds-training/sep20/shared/week2/ERR1755082.test.sort.bam' 

parser = argparse.ArgumentParser(prog='This is a BAM to BED conversion')
parser.add_argument('--input', '-i', dest='bamfilepath', help='Input bam file')
parser.add_argument('--output', '-o', dest='outputfilepath', help='Output bed file')
args = parser.parse_args()

#/Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/obds_training/Python/ERR1755082.test.sort.bam' 
bamfile = pysam.AlignmentFile("args.bamfilepath", "rb")
bamfile_iter = bamfile.fetch() 
outputfilepath = "/Users/Sophie/Documents/UniversityofOxford/DPhil/OBDS/obds_training/Python/outputfile.txt"
with open(outputfilepath, 'w') as outputfile:
    for alignmnet in bamfile_iter:
        if alignmnet.is_paired:
            #/n = make a new line 
            outputfile.write(f'''{alignment_reference_name}\t{alignment_reference_start}\t{query_alignment_end},t{query_name}\t{alignment.mapping_quality}\t.\n''')
        
            # print(f'''{alignment_reference_name}\t{alignment_reference_start}\t{query_alignment_end},t{query_name}\t{alignment.mapping_quality}\t.''')

        
       # for column in bamfile(contig, start, end, truncate = truncate):
       #      if truncate:
       #          self.assertGreaterEqual( column.pos, start )
       #          self.assertLess( column.pos, end )
            
            
            # samfile.close()

            
            
            
  # with pysam.sort("-o", "output.bam", "ERR1755082.bam") :            