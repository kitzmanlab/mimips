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

    opts.add_argument('--qileCvgs', default='0.10,0.25,0.50,0.75,0.90', dest='qileCvgs',
        help='comma sepd quantiles over which to check cvg, eg 0.25,0.50,0.75')

    opts.add_argument('--cvgThresholds', default=None, dest='cvgThresholds')    

    opts.add_argument('--coloutPrefix',default="",dest='coloutPrefix')

    opts.add_argument('--outKey',default=None,dest='outKey')

    o = opts.parse_args()

    ######################################################

    lThresholds = [] if o.cvgThresholds is None else [ int(x) for x in o.cvgThresholds.split(',') ]

    lQiles = [] if o.qileCvgs is None else [  float(x) for x in o.qileCvgs.split(',') ]

    tblIn = pd.read_csv(o.inKey,sep='\t')

    for thresh in lThresholds:
        tblIn['%scvg_gte_%d'%( o.coloutPrefix, thresh )] = 0

    for qile in lQiles:
        tblIn['%scvg_qile_%02d'%( o.coloutPrefix,int(100*qile) ) ] = 0

    for i,r in tblIn.iterrows():

        tblCvg = pd.read_csv( r[o.colHisto], sep='\t' )

        for thresh in lThresholds:
            fracOverThresh = tblCvg.ix[ tblCvg['depth']>=thresh, 'percent' ].sum()

            tblIn.loc[i,'%scvg_gte_%d'%( o.coloutPrefix, thresh )] = fracOverThresh

        for qile in lQiles:
            print(qile)
            liAbvQile = tblCvg['cumulative_percent'] >= qile*100
            if not liAbvQile.any():
                tblIn.loc[i, '%scvg_qile_%02d'%( o.coloutPrefix, int(100*qile) ) ] = 0
            else:
                qilepos=liAbvQile[liAbvQile].index[-1]
                tblIn.loc[i, '%scvg_qile_%02d'%( o.coloutPrefix, int(100*qile) ) ] = int(tblCvg.loc[qilepos,'depth'])


    tblIn.to_csv(o.outKey,sep='\t',index=False)

if __name__ == '__main__':                
    main()