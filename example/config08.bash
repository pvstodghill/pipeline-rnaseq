#! /bin/bash

run_deseq2 -f 1.0 -c 2.0 -t example WT:WT1 WT:WT2 WT:WT3 5255:52551 5255:52552 5255:52553

run_deseq2 -f -1 -c 1.5 -t unfiltered WT:WT1 WT:WT2 WT:WT3 5255:52551 5255:52552 5255:52553
