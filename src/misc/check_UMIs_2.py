#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This code checks the total number of detected UMIs and CBCs in a tagged bam file.

The expected format of the tags is a sfollows:
XM: UMI
XC: cell barcode (CBC), a concatenation of the three splitseq barcodes
gn: gene symbol, as attached by the Dropseqtools TagReadWithInterval and TagReadWithGeneFunction

@author: rebekka
"""

#%% libraries
from __future__ import division 
import sys, re, os
import pysam
import numpy as np

#%% Functions
def hamming(s1, s2):
    """Calculate the Hamming distance between two strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

#%% Read in BAM file
#print(sys.argv[1])
# note: if you provide "-" as an argument, this can read from stdin
    
infile = pysam.AlignmentFile(sys.argv[1], 'r', check_sq=False, check_header=False)
output_file = 'read_vs_umi.txt'

if not os.path.isfile(output_file):
    print('Creating new ouput file')
    with open(output_file, 'w') as f:
        print('No.reads\tNo.UMIs\tNo.CBCs', file = f)
      
n=0
n1=0


UMIs = dict()
CBCs = dict()

for entry in infile.fetch(until_eof=True):
    if not entry.has_tag('gn') or entry.mapq < 3:
        n1+=1
        continue
    n+=1
    
    xm = entry.get_tag('XM')
    gn = entry.get_tag('gn')
    xc = entry.get_tag('XC')

    if xc in CBCs:
        CBCs[xc]+=1
    else:
        CBCs.update({xc:1})

    combo = xc+xm+gn

    if xm in UMIs:
        UMIs[combo]+=1
    else:
        UMIs.update({combo:1})

with open(output_file, 'a') as f:
    print('%d\t%d\t%d' % (n,len(UMIs.keys()),len(CBCs.keys())), file = f)

infile.close()

print('skipped %d lines' % n1)
