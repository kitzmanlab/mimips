import sys
import argparse
from collections import defaultdict, OrderedDict

import os.path

import pysam

import pandas as pd

def main():    
    opts = argparse.ArgumentParser()
    
    opts.add_argument('--inKey',default=None,dest='inKey')
    
    opts.add_argument('--colPerBase',default='cvg_perbase',dest='colPerBase')

    opts.add_argument('--cvgThresholds', default=None, dest='cvgThresholds')    

    opts.add_argument('--outTbl',default=None,dest='outTbl')

    o = opts.parse_args()

    ######################################################

    lThresholds = [] if o.cvgThresholds is None else [ int(x) for x in o.cvgThresholds.split(',') ]
    assert len(lThresholds)>=1

    tblIn = pd.read_csv(o.inKey,sep='\t')


    tblOut = pd.DataFrame( OrderedDict(
        [ (c,[]) for c in ['tgtchrom','tgtstart','tgtend','pos']+['pct_samps_gte_%d'%(thresh) for thresh in lThresholds]
             ] ) )


    ctr=0

    for i,r in tblIn.iterrows():

        tblCur = pd.read_csv( r[o.colPerBase], sep='\t' )

        sys.stderr.write('gather per-bp coverage from file %d/%d\n'%( ctr, tblIn.shape[0] ))
        sys.stderr.flush()

        if ctr==0:
            # get positions from this first file
            tblOut['tgtchrom']=tblCur['tgtchrom']
            tblOut['tgtstart']=tblCur['tgtstart']
            tblOut['tgtend']=tblCur['tgtend']
            tblOut['pos']=tblCur['pos']

            for c in ['tgtstart','tgtend','pos']:
                tblOut[c]=tblOut[c].astype('int')

        for thresh in lThresholds:
            if ctr==0:
                tblOut['pct_samps_gte_%d'%thresh] = (tblCur['cvg'] >= thresh).astype('int')
            else:
                tblOut['pct_samps_gte_%d'%thresh] += (tblCur['cvg'] >= thresh).astype('int')

        ctr+=1

    for thresh in lThresholds:
        tblOut['pct_samps_gte_%d'%thresh] /= float(ctr)
        tblOut['pct_samps_gte_%d'%thresh] *= 100

    tblOut.to_csv( o.outTbl, sep='\t', index=False )        



if __name__ == '__main__':                
    main()