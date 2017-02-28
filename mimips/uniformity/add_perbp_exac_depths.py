import sys
import argparse
from collections import defaultdict, OrderedDict

import os.path

import pysam

import numpy as np

import pandas as pd

from cStringIO import StringIO

def main():    
    opts = argparse.ArgumentParser()
    
    opts.add_argument('--inPerbpMIPS',default=None,dest='inPerbpMIPS')
    opts.add_argument('--inPerbpExac',default=None,dest='inPerbpExac')

    opts.add_argument('--outPerbpJoined',default=None,dest='outPerbpJoined')
    
    o = opts.parse_args()

    ######################################################

    tblIn = pd.read_table(o.inPerbpMIPS)
    tblIn['tgtchrom']=tblIn['tgtchrom'].astype('str')

    tblOut = tblIn.set_index( ['tgtchrom','pos'] )
    
    for c in ['exac_mean','exac_median']+['exac_pctat_%d'%cvg for cvg in (1,5,10,15,20,25,30,50,100)]:
        tblOut[c]=-1

    tbxExac = pysam.TabixFile( o.inPerbpExac )

    for chrom in tblIn['tgtchrom'].unique():
        print chrom

        try:
            curChromExac = tbxExac.fetch( chrom )
        except:
            continue

        sio = StringIO('\n'.join(curChromExac))
        tblChromExac = pd.read_table( sio, header=None )
        tblChromExac.columns = ['chrom','pos','exac_mean','exac_median']+\
                    ['exac_pctat_%d'%cvg for cvg in (1,5,10,15,20,25,30,50,100)]
        del tblChromExac['chrom']

        for c in ['exac_pctat_%d'%cvg for cvg in (1,5,10,15,20,25,30,50,100)]:
            tblChromExac[c]*=100.


        tblInCurCurom = tblIn.query("tgtchrom=='%s'"%chrom)

        tblInCurChromJoin = pd.merge( tblInCurCurom, tblChromExac, left_on='pos', right_on='pos', how='left' )
    
        for c in ['exac_mean','exac_median']+['exac_pctat_%d'%cvg for cvg in (1,5,10,15,20,25,30,50,100)]:
            tblOut.loc[  zip(tblInCurChromJoin['tgtchrom'],tblInCurChromJoin['pos']), c ] = np.array(tblInCurChromJoin[c])


    tblOut.to_csv( o.outPerbpJoined, sep='\t', index=True, na_rep='NA' )


if __name__=='__main__':
    main()