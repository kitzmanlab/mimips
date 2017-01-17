from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
import sys
import argparse
from collections import defaultdict

import os.path

import pysam

import numpy as np

import pandas as pd

from itertools import groupby

from mip_pipe_common import *

import math
def f(x,c,n):
  return old_div(c,x) - 1 + math.exp(old_div(-n,x))


def estLibSize( ntot, ndup )    :

    if ntot == 0:
        return -1

    duprate = old_div(ndup,float(ntot))

    if ntot > 1000 and duprate > 0: 

        c=int( (1-duprate )*ntot)

        m=1.
        M=100.

        while ( f(M*c, c, ntot) >= 0 ):
          M *= 10.0

        for i in range(40):
            r = old_div((m+M),2.0);
            u = f( r * c, c, ntot );
            if ( u == 0 ) : break
            elif ( u > 0 ) : m = r
            elif ( u < 0 ) : M = r

        OUT=c*(m+M)/2.

        return OUT
    else:
        return -1 

def main():

    opts = argparse.ArgumentParser()

    # opts.add_argument('--pass1',default=False, action='store_true', dest='pass1',
    #                    description='pass 1 : reads bam ')

    opts.add_argument('--libname',default='not_specified',dest='libname')

    opts.add_argument('--mergedReads',default=False,action='store_true',dest='mergedReads')

    opts.add_argument('--inBam',default=None,dest='inBam')
    opts.add_argument('--outBam',default=None,dest='outBam')

    opts.add_argument('--suppressDupMIDs',default=False,action='store_true',dest='suppressDupMIDs')

    opts.add_argument('--outSummary',default='/dev/null',dest='outSummary')
    
    o = opts.parse_args()

    ######################################################

    # pass 1 - run on pos sorted bam
    #  cache readnames that are deemed MID-dups
    # pass 2 - reiterate pos sorted bam and flag those reads and their mates (wherever they may be)

    mergedReads = o.mergedReads

    bamIn = pysam.Samfile( o.inBam, 'rb' )

    sReadidsAreDup = set()

    nProcessedPass1=0
    nPass1AreDup=0
    nMidsAreDup=0

    for chrompos, lReads in groupby( bamIn, lambda r:(r.reference_id, r.reference_start ) ):
        
        if mergedReads:
            lr1 = sorted( lReads, key=lambda r: dict(r.tags)['YT'] )    
        else:
            lr1 = sorted( [ r for r in lReads if r.is_read1 ], key=lambda r: dict(r.tags)['YT'] )

        for _,lr1_bytag in groupby(lr1, lambda r: dict(r.tags)['YT']):

            lr1_bytag = list(lr1_bytag)

            # more than one read at this position with this MID?
            if len(lr1_bytag) > 1:
                lr1_bytag = sorted( lr1_bytag, key=lambda r: -np.mean(r.query_alignment_qualities) )
                sReadidsAreDup.update( [ r1.query_name for r1 in lr1_bytag[1:] ] )
                nPass1AreDup+=len(lr1_bytag)-1
                nMidsAreDup+=1

            nProcessedPass1 += len(lr1_bytag)

            if (old_div((nProcessedPass1-len(lr1_bytag)),10000)) < (old_div(nProcessedPass1,10000)):
                print('processed %d read 1s'%(nProcessedPass1)) 
                print('identified %d mid duplicates accounting for %d reads'%(nMidsAreDup, nPass1AreDup))
                print('total dup read rate = %.2f%%'%( 100.*nPass1AreDup/float(1+nProcessedPass1) )) 

    print('processed %d read 1s'%(nProcessedPass1)) 
    print('identified %d mid duplicates accounting for %d reads'%(len(sReadidsAreDup), nPass1AreDup))
    print('total dup read rate = %.2f%%'%( 100.*nPass1AreDup/float(1+nProcessedPass1) )) 
    totalLibSize = estLibSize( nProcessedPass1, nPass1AreDup )
    print('total libsize = %.2e'%( totalLibSize ))

    bamIn.reset()

    bamOut = pysam.Samfile( o.outBam, 'wb', template=bamIn )

    nFlaggedRead1s = 0

    if not mergedReads:
        if not o.suppressDupMIDs:
            for r in bamIn:
                if r.query_name in sReadidsAreDup:
                    if r.is_read1:
                        nFlaggedRead1s+=1
                    r.flag|=0x400
                bamOut.write(r)    
        else:
            for r in bamIn:
                if r.query_name in sReadidsAreDup:
                    if r.is_read1:
                        nFlaggedRead1s+=1
                else:
                    bamOut.write(r)    
    else:
        if not o.suppressDupMIDs:
            for r in bamIn:
                if r.query_name in sReadidsAreDup:
                    nFlaggedRead1s+=1
                    r.flag|=0x400
                bamOut.write(r)    
        else:
            for r in bamIn:
                if r.query_name in sReadidsAreDup:
                    nFlaggedRead1s+=1
                else:
                    bamOut.write(r)    


    filOutSummary = open(o.outSummary,'w')
    filOutSummary.write( tabnl( ['libname','duprate','n_read1_checked_for_dup_mids','n_dup_mids_found','est_lib_size'] ))
    filOutSummary.write( tabnl( [o.libname, 
                                 '%.3f'%(old_div(float(len(sReadidsAreDup)),(1+nProcessedPass1))),
                                 nProcessedPass1,
                                 len(sReadidsAreDup),
                                 '%.3e'%(totalLibSize) ] ) )

if __name__ == '__main__':                
    main()