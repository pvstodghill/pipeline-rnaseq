#! /bin/bash

# these are the default cutoffs
run_deseq2 -f 0.05 -c 2.0 -t default WT:WT1 WT:WT2 WT:WT3 5255:52551 5255:52552 5255:52553

# these cuttoffs are tight
run_deseq2 -f 0.01 -c 2.0 -t publishable WT:WT1 WT:WT2 WT:WT3 5255:52551 5255:52552 5255:52553

# don't apply FDR cutoff at all; lower fold-change threshold to 1.5
run_deseq2 -f -1 -c 1.5 -t unfiltered WT:WT1 WT:WT2 WT:WT3 5255:52551 5255:52552 5255:52553
