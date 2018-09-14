#! /usr/bin/env python

from __future__ import print_function
import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import numpy as np
import numpy.random as rand
import re

import pandas as pd
import pybedtools as pbt

from ..mip_pipe_common import *

def main():    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inPerProbeUnifTbl', dest='inPerProbeUnifTbl')
    opts.add_argument('--libname', dest='libname')
    opts.add_argument('--bedXY', dest='bedXY',default=None,help='bed file containing specific X and Y regions to check')
    opts.add_argument('--outSexChromTbl', dest='outSexChromTbl')

    o = opts.parse_args()

    if o.bedXY is not None:
        btxy = pbt.BedTool( o.bedXY )
    else:
        btxy = None

    tblInMips = pd.read_table(o.inPerProbeUnifTbl)
    
    sumhits = tblInMips['probe_hits_OK'].sum()

    if btxy is not None:
        tblInMips_prefilt = tblInMips
        tblInMips = bed_intersect_dataframe( 
            btxy, 
            tblInMips_prefilt)

    if sumhits > 0:
        tblInMipsXY = tblInMips.loc[ (tblInMips.chrom=='X')|(tblInMips.chrom=='Y') ]
        sumhitsXY = tblInMipsXY[ ['chrom','probe_hits_OK'] ].groupby('chrom').agg(sum)
        frachitsXY =  sumhitsXY / float(sumhits)
        filOut=open(o.outSexChromTbl,'w')
        filOut.write('libname\tchrX_frac\tchrY_frac\n')
        filOut.write('{}\t{:.3e}\t{:.3e}\n'.format( 
            o.libname, 
            frachitsXY.loc['X','probe_hits_OK'], 
            frachitsXY.loc['Y','probe_hits_OK'] ))
        filOut.close()
    else:
        filOut=open(o.outSexChromTbl,'w')
        filOut.write('libname\tchrX_frac\tchrY_frac\n')
        filOut.write('{}\t0\t0\n')
        filOut.close()


if __name__ == '__main__':
    main()