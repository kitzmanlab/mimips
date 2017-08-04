#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
# #from builtins import str  # this breaks all sorts of other scripts which test type(x)==str    # this breaks all sorts of other scripts which test type(x)==str
from past.utils import old_div
import sys

import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import numpy as np
import numpy.random as rand
import re

import pandas as pd

import pysam

from jkutils import *

if __name__ == '__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--mipTableIn', action='append', dest='mipTableIn')

    opts.add_argument('--inDesignName', action='append', dest='inDesignName')

    opts.add_argument('--mipTableOut', default=None, dest='mipTableOut')

    o = opts.parse_args()

    """

    #specify --mipTableIn file --inDesignName name for each design to be concatenated

    join_multi_design_tables.py --mipTableIn design1.pairs.txt --inDesignName design1 --mipTableIn design2.pairs.txt --inDesignName design2 [ ... ] --mipTableOut

    """

    idesi=0
    mDesiToOuterCorng={}
    for fn in o.mipTableIn:
        desi=o.inDesignName[idesi]
        mDesiToOuterCorng[desi]=set()
        tblCur = pd.read_csv(fn,sep='\t')
        for _,r in tblCur.iterrows():
            mDesiToOuterCorng[desi].add( (r.chrom,r.start,r.end) )
        idesi+=1

    sOuterAlreadyRecorded=set()

    lTables = [ ]
    idesi=0
    idxBase=0
    for fn in o.mipTableIn:
        lTables.append( pd.read_csv(fn,sep='\t') )
            
        desi=o.inDesignName[idesi]

        lDesiMultiNames=[]

        for _,r in lTables[-1].iterrows():
            lDesiMultiNames.append('_'.join( [ o.inDesignName[jdesi] for jdesi in range(len(o.inDesignName)) 
                                      if (r.chrom,r.start,r.end) in mDesiToOuterCorng[o.inDesignName[jdesi]] ] ))

        lTables[-1]['design']=lDesiMultiNames

        lKeepEntry = [ (r.chrom,r.start,r.end) not in sOuterAlreadyRecorded 
                       for _,r in lTables[-1].iterrows() ]

        lTables[-1] = lTables[-1].loc[lKeepEntry]

        sOuterAlreadyRecorded = sOuterAlreadyRecorded.union(
            [ (r.chrom,r.start,r.end) for _,r in lTables[-1].iterrows() ] )

        lTables[-1]['mip_index']=lTables[-1]['mip_index']+idxBase
        idesi+=1
        idxBase=lTables[-1]['mip_index'].max()+1

    mipTblOut=pd.concat(lTables)

    mipTblOut.to_csv(o.mipTableOut,sep='\t',index=False)





