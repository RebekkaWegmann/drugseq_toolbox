#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is to process DRUG-seq reads after tagging the bam file with well and molecular barcodes.

Input: BAM file with reads that have been filtered, polyA and SMART adapter trimmed, and tagged with the following:
    - XW: well barcode
    - XM: molecular barcode (UMI)

Output:
    - Filtered bam file containing only reads with expected barcodes. Barcodes that are within 1 hamming distance of an expected
    barcode sequence are corrected. The plate position is added to the well barcode tag.

Copyright: Rebekka Wegmann, Snijderlab, ETH Zurich, 2019

# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
    
"""

#%% libraries
import sys, os, timeit, h5py
import pysam
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf') #prevents matplotlib from trying to open a figure window, which will fail on systems that do not have graphics
from matplotlib import pyplot as plt
from argparse import ArgumentParser

#%% Functions
def hamming(s1, s2):
    """Calculate the Hamming distance between two strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def well2ind(well):
    """Convert well positions to (zero-based) indices"""
    d = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7}
    row = d[well[0]]
    col = int(well[1:])-1
    return [row,col]

def make_plate_overview(bc_list, bc_count):
    """Convert the list of barcode counts to counts per well for plotting"""
    out = np.zeros([8,12])
    bc2ind = {bc_list.Barcode.values[x]:well2ind(bc_list.WellPosition.values[x]) for x in range(len(bc_list.Barcode.values))}
    for bc in bc2ind.keys():
        out[bc2ind[bc][0],bc2ind[bc][1]] = bc_count[bc]
    return out

#%%Setup input parser
parser = ArgumentParser()
parser.add_argument("-i" "--input_bam", action="store", dest="input_bam", default="-", help="Specify the input bam file. Defaults to stdin.")
parser.add_argument("-o" "--output_bam", action="store", dest="output_bam", default="-", help="Specify the output bam file. Defaults to stdout.")
parser.add_argument("-d" "--out_dir", action="store", dest="out_dir", default=".", help="Directory to store logfiles and output plots. Defaults to the current directory.")
parser.add_argument("-b" "--bc_dir", action="store", dest="bc_dir", default=".", help="Directory where the expected barcode files are stored. Defaults to the directory this script is in.")
parser.add_argument("--debug_flag",action="store_true",help="Turn on debug flag. This will produce some additional output which might be helpful.")
parser.add_argument("--store_discarded",action="store_true",help="Store names of discarded reads?")

#%% Parse input
args = parser.parse_args()

debug_flag = args.debug_flag
store_discarded = args.store_discarded
input_bam = args.input_bam #use "-" for stdin, set flag to rb
output_bam = args.output_bam
out_dir = args.out_dir
bc_dir = args.bc_dir

if bc_dir==".":
    bc_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
    
#%% Write parameters to logfile
print('Drugseq barcode filtering log\n---------------------------------------\nParameters:', file = open(os.path.join(out_dir,'barcode_filtering_log.txt'), 'w'))
print('Input bam: %s' % input_bam, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Ouput bam: %s' % output_bam, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Output directory: %s' % out_dir, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Barcode directory: %s' % bc_dir, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))

if store_discarded:
    if os.path.isfile(os.path.join(out_dir,'discarded_reads.txt')):
        os.remove(os.path.join(out_dir,'discarded_reads.txt'))
        print('Old version of discarded reads.txt deleted.')

#%% Start timing
t_start = timeit.default_timer()

#%% Get expected barcodes. 
expected_bc = pd.read_csv(os.path.join(bc_dir,"expected_barcodes.csv"))

#%% Read in BAM file
infile =pysam.AlignmentFile(input_bam, 'rb', check_sq=False)
outfile = pysam.AlignmentFile(output_bam, 'wb', template=infile)

#%% for debugging only
if debug_flag:
    for entry in itertools.islice(infile, 10 ):
        print(entry.query_name)
        print(entry.get_forward_sequence())
        print(entry.get_tag('XW'))
        print(entry.get_tag('XM'))
        
#%% Check and correct barcodes

n_entries=0
n_correct_xw=0
n_edited_bc=0

n_bc  = {x:0 for x in expected_bc.Barcode.values}  
all_bcs =  [None]*50000


for entry in infile.fetch(until_eof=True):
    n_entries+=1
    keep=True
    xw = entry.get_tag('XW')

    if n_entries<50000:
    	all_bcs[n_entries]=xw
    
    if xw in expected_bc.Barcode.values:
        n_correct_xw+=1
        n_bc[xw]+=1

    else:
        d = [hamming(xw,bc) for bc in expected_bc.Barcode.values]
        idx = [i for i,e in enumerate(d) if e==1]
        if len(idx)==1:
            n_edited_bc+=1
            xw = expected_bc.Barcode.values[idx[0]]
            entry.set_tag('XW',xw)
            n_bc[xw]+=1
        else:
            keep=False
           
    if keep:
        well = expected_bc.WellPosition.values[expected_bc.Barcode.values==xw][0]
        xw = well+'_'+xw
        entry.set_tag('XW',xw)
        outfile.write(entry)
    else:
        if store_discarded:
            print(entry.query_name, file = open(os.path.join(out_dir,'discarded_reads.txt'),'a'))

print('Read %d entries' % n_entries, file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Found %d [%.2f%%] correct barcodes' % (n_correct_xw, n_correct_xw/float(n_entries)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Corrected %d [%.2f%%] barcodes' % (n_edited_bc , n_edited_bc/float(n_entries)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))
print('Retained %d [%.2f%%] reads after correction and filtering' % (n_edited_bc+n_correct_xw, (n_edited_bc+n_correct_xw)/float(n_entries)*100), file = open(os.path.join(out_dir,'barcode_filtering_log.txt'),'a'))

infile.close()
outfile.close()

t_stop = timeit.default_timer()
t_elapsed = t_stop-t_start
print("Barcode filtering finished. Elapsed time: %ds" % t_elapsed)
t_start=t_stop

#%% Create summary plots
    
bc_matrix = make_plate_overview(expected_bc, n_bc)

# this part is super slow if you have a large number of cells, maybe omit
all_bc_counts = {i:all_bcs.count(i) for i in list(set(all_bcs))}
all_bc_cumsum = np.cumsum(sorted(list(all_bc_counts.values()), reverse=True))/float(50000)

fig, axes = plt.subplots(1,2)
p1 = axes[0].imshow(np.log10(bc_matrix+1));  axes[0].set_title('Number of reads per well barcode')
clb = fig.colorbar(p1, ax=axes[0]); clb.set_label('No. BCs, log10')
p2 = axes[1].plot(all_bc_cumsum); axes[1].set_title('Cumulative fraction of reads per barcode')

fig.set_size_inches(10,4)
fig.savefig(os.path.join(out_dir,'barcode_filtering_summary.pdf'),bbox_inches='tight')

#%% Stop timing
t_stop = timeit.default_timer()
t_elapsed = t_stop-t_start

#%% Save the QC output
f = h5py.File(os.path.join(out_dir,'splitseq_filtering_QC_data.hdf5'), 'w')

f.create_dataset('plate_overview', data = bc_matrix)
f.create_dataset('reads_per_BC', data=list(all_bc_counts.values()))
f.create_dataset('labels_reads_per_BC', data=np.string_(list(all_bc_counts.keys())))

f.close()

#%% print info to stdout
print("Generating summary finished. Elapsed time: %ds" % t_elapsed)
