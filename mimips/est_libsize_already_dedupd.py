from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
import sys
import argparse
from collections import defaultdict

import os.path

import pysam

import pandas as pd

from jkutils import *

from itertools import groupby


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

def main() :   
    opts = argparse.ArgumentParser()

    opts.add_argument('--libname',default='not_specified',dest='libname')

    opts.add_argument('--inPreConsensusBam',default=None,dest='inPreConsensusBam')
    opts.add_argument('--inPostConsensusBam',default=None,dest='inPostConsensusBam')

    opts.add_argument('--outSummary',default='/dev/null',dest='outSummary')
    
    o = opts.parse_args()

    ######################################################

    lr1 = pysam.flagstat(o.inPreConsensusBam).split('\n')[6].strip()
    assert lr1.endswith('read1')
    nPreCns = int(lr1.split(' ')[0])

    lr1 = pysam.flagstat(o.inPostConsensusBam).split('\n')[6].strip()
    assert lr1.endswith('read1')
    nPostCns = int(lr1.split(' ')[0])

    print('processed %d read 1s'%(nPreCns)) 
    print('total dup read rate = %.2f%%'%( 100.*(nPreCns-nPostCns)/(1+nPreCns)) )
    totalLibSize = estLibSize( nPreCns, nPreCns-nPostCns )
    print('total libsize = %.2e'%( totalLibSize ))

    filOutSummary = open(o.outSummary,'w')
    filOutSummary.write( tabnl( ['libname','duprate','n_read1_checked_for_dup_mids','n_dup_mids_found','est_lib_size'] ))
    filOutSummary.write( tabnl( [o.libname, 
                                 '%.3f'%( (1.0 - (nPostCns/(1+nPreCns)) )),
                                 nPreCns,
                                 nPreCns-nPostCns,
                                 '%.3e'%(totalLibSize) ] ) )

if __name__ == '__main__':                
    main()