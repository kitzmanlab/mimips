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

import pysam

from ..mip_pipe_common import *

def main():    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--inMipPairTbl', dest='inMipPairTbl')
    opts.add_argument('--inBam', dest='inBam')
    opts.add_argument('--outMipPairTbl', dest='outMipPairTbl')
    opts.add_argument('--prefixChar', default='Z', dest='prefixChar')

    o = opts.parse_args()

    tblInMips = pd.read_csv(o.inMipPairTbl,sep='\t')
    tblInMips = tblInMips.set_index('mip_index')

    tblInMips['ext_arm_only_hits'] = 0
    tblInMips['lig_arm_only_hits'] = 0

    tblInMips['ext_arm_mispaired_hits'] = 0
    tblInMips['lig_arm_mispaired_hits'] = 0

    tblInMips['probe_hits_acap'] = 0

    tblInMips['probe_hits_OK'] = 0

    ctr=0

    mProbeIdHits=defaultdict( lambda:{'probe_hits_OK':0,'probe_hits_acap':0} )
    mArmIdHits=defaultdict( lambda:{'ext_arm_only_hits':0,'lig_arm_only_hits':0,
                                   'ext_arm_mispaired_hits':0, 'lig_arm_mispaired_hits':0} )


    tagp,tage,tagl,tagc='%sP'%o.prefixChar, '%sE'%o.prefixChar, '%sL'%o.prefixChar, '%sC'%o.prefixChar

    bam = pysam.Samfile( o.inBam, 'rb' )
    try:
        for rd in bam:

            if ctr%10000 == 0: print('%d...'%ctr)

            ctr+=1

            tags = dict(rd.tags)

            yp = tags[tagp]
            ye = tags[tage]
            yl = tags[tagl]
            yc = tags[tagc]

            if yp >= 0:
                if yc != 0:
                    # tblInMips.ix[ yp, 'probe_hits_acap' ]  += 1  # ...much slower...
                    mProbeIdHits[yp]['probe_hits_acap']+=1
                else:
                    # tblInMips.ix[ yp, 'probe_hits_OK' ]  += 1
                    mProbeIdHits[yp]['probe_hits_OK']+=1
            else:
                if ye>=0 and yl<0:
                    # tblInMips.ix[ ye, 'ext_arm_only_hits' ] += 1
                    mArmIdHits[ye]['ext_arm_only_hits']+=1
                elif yl>=0 and ye<0:
                    # tblInMips.ix[ yl, 'lig_arm_only_hits' ] += 1
                    mArmIdHits[yl]['lig_arm_only_hits']+=1
                elif ye>=0 and yl>=0:
                    # tblInMips.ix[ ye, 'ext_arm_mispaired_hits' ] += 1
                    # tblInMips.ix[ yl, 'lig_arm_mispaired_hits' ] += 1
                    mArmIdHits[ye]['ext_arm_mispaired_hits']+=1
                    mArmIdHits[yl]['lig_arm_mispaired_hits']+=1


    except:
        print(sys.exc_info())
    
    tblInMips['probe_hits_OK'] = [ mProbeIdHits[yp]['probe_hits_OK'] for yp in tblInMips.index ]
    tblInMips['probe_hits_acap'] = [ mProbeIdHits[yp]['probe_hits_acap'] for yp in tblInMips.index ]

    tblInMips['ext_arm_only_hits'] = [ mArmIdHits[yp]['ext_arm_only_hits'] for yp in tblInMips.index ]
    tblInMips['lig_arm_only_hits'] = [ mArmIdHits[yp]['lig_arm_only_hits'] for yp in tblInMips.index ]
    tblInMips['ext_arm_mispaired_hits'] = [ mArmIdHits[yp]['ext_arm_mispaired_hits'] for yp in tblInMips.index ]
    tblInMips['lig_arm_mispaired_hits'] = [ mArmIdHits[yp]['lig_arm_mispaired_hits'] for yp in tblInMips.index ]

    tblInMips.reset_index(drop=False).to_csv( o.outMipPairTbl,sep='\t', index=False )

if __name__ == '__main__':
    main()