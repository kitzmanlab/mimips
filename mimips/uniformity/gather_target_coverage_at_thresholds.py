import sys
import argparse
from collections import defaultdict

import os.path

import pysam

import pandas as pd

import pybedtools as pbt

from ..mip_pipe_common import *

def main():    
    opts = argparse.ArgumentParser()
    
    opts.add_argument('--inKey',default=None,dest='inKey')
    
    opts.add_argument('--colHisto',dest='colHisto')

    opts.add_argument('--cvgThresholds', default=None, dest='cvgThresholds')    

    opts.add_argument('--coloutPrefix',default="",dest='coloutPrefix')

    opts.add_argument('--outKey',default=None,dest='outKey')

    o = opts.parse_args()

    ######################################################

    lThresholds = [] if o.cvgThresholds is None else [ int(x) for x in o.cvgThresholds.split(',') ]

    tblIn = pd.read_csv(o.inKey,sep='\t')

    for thresh in lThresholds:
        tblIn['%scvg_gte_%d'%( o.coloutPrefix, thresh )] = 0

    for i,r in tblIn.iterrows():

        tblCvg = pd.read_csv( r[o.colHisto], sep='\t' )

        for thresh in lThresholds:
            fracOverThresh = tblCvg.ix[ tblCvg['depth']>=thresh, 'percent' ].sum()

            tblIn.ix[i,'%scvg_gte_%d'%( o.coloutPrefix, thresh )] = fracOverThresh

    tblIn.to_csv(o.outKey,sep='\t',index=False)

if __name__ == '__main__':                
    main()