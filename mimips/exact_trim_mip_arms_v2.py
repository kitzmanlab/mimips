from builtins import range
import sys
import argparse
from collections import defaultdict
from itertools import groupby
from math import copysign

import os.path

import pysam
import pandas as pd

from mip_pipe_common import *
from mip_pipe_core import *

def clipReadToGapfill( gapFillStart, gapFillEnd, rd ):
    
    """
    gapFillStart = 0-based coordinate of gap fill start
    gapFillEnd = 0-based coordinate of gap fill end

    returns:
        read, modified
    """

    cigTups = rd.cigartuples

    rdq = rd.query_qualities
    rds = rd.query_sequence

    # easy case: existing cigar str is just a single match block
    if len(cigTups) == 1 and cigTups[0][0]==0:

        # rd.ref_end : 0, non-inclusive
        cozrefEndClipToGapFill = min( rd.reference_end-1, gapFillEnd )

        rnginclQueryNew = \
            ( max( rd.reference_start, gapFillStart ) - rd.reference_start,
              min( rd.reference_end-1, gapFillEnd ) - rd.reference_start )

        rd.cigartuples = ( (0,rnginclQueryNew[1]-rnginclQueryNew[0]+1), )

        rd.reference_start = rd.reference_start + rnginclQueryNew[0]
        rd.query_sequence = rds[ rnginclQueryNew[0]:rnginclQueryNew[1]+1 ]
        rd.query_qualities = rdq[ rnginclQueryNew[0]:rnginclQueryNew[1]+1 ]

        return rd


    # if this fails, modify and then rebuild the cigar string.
    nbpCig = sum( [ ci[1] for ci in cigTups ] )

    lNongapRdRef = -999 * np.ones( (nbpCig,2), 'int32' )
    cozgNext = 0
    for ci in rd.cigartuples:
        if ci[0]==0:
            lNongapRdRef[ cozgNext:cozgNext+ci[1], : ] = 1 
        elif ci[0]==1: # insertion to ref
            lNongapRdRef[ cozgNext:cozgNext+ci[1], 0 ] = 1
            lNongapRdRef[ cozgNext:cozgNext+ci[1], 1 ] = 0
        elif ci[0]==2: # deletion from ref
            lNongapRdRef[ cozgNext:cozgNext+ci[1], 0 ] = 0
            lNongapRdRef[ cozgNext:cozgNext+ci[1], 1 ] = 1
        elif ci[0]==4: # masked query
            lNongapRdRef[ cozgNext:cozgNext+ci[1], 0 ] = -1
            lNongapRdRef[ cozgNext:cozgNext+ci[1], 1 ] = 0
        else: assert 1==0, 'other cigar codes not supported.'
        cozgNext += ci[1]
    
    assert (lNongapRdRef==-999).sum()==0

    # find the FIRST aligned position AT OR INSIDE the gap fill
    ixrefGfstart = max( rd.reference_start, gapFillStart ) - rd.reference_start

    # find the LAST aligned position AT OR INSIDE the gap fill
    ixrefGfend = min(rd.reference_end - 1, gapFillEnd) - rd.reference_start

    #        0  1  2  3  4  5  6  7  8  9 10 11
    #
    #  rd    1  1  1  1  0  1  1  1  1  1  1  1  
    # csrd   1  2  3  4  4  5  6  7  8  9 10 11
    #  ref   1  1  0  1  1  1  1  1  1  1  1  1
    # csref  1  2  2  3  4  5  6  7  8  9 10 11
    # refcz 
    #  100+  0  1  1  2  3  4  5  6  7  8  9 10
    #                       *              *  
    #                   
    #    refstart=100  refend=111   
    #    gfstart=104     ixrefGfstart = 4    ixgapalnGfstart = 5
    #    gfend=109       ixrefGfend = 9      ixgapalnGfend = 10
    #    
    #  

    csNongapRdRef = (lNongapRdRef != 0).cumsum(0)

    ixgapalnGfstart = csNongapRdRef[:,1].searchsorted( ixrefGfstart + 1,  'left' )
    ixgapalnGfend = csNongapRdRef[:,1].searchsorted( ixrefGfend + 1,  'left' )

    # make sure we do not start @ deletion
    while ixgapalnGfstart < nbpCig and lNongapRdRef[ ixgapalnGfstart, 0 ] == 0:
        ixgapalnGfstart += 1

    if ixgapalnGfstart == nbpCig or ixgapalnGfstart >= ixgapalnGfend:
        return None

    # make sure we do not end @ deletion
    while ixgapalnGfend >= 0 and lNongapRdRef[ ixgapalnGfend, 0 ] == 0:
        ixgapalnGfend -= 1

    if ixgapalnGfend < 0:
        return None

    ################ 

    cozrefNewStart = rd.reference_start + csNongapRdRef[ ixgapalnGfstart , 1 ] - 1

    rnginclQueryNew = ( csNongapRdRef[ ixgapalnGfstart, 0 ] - 1,
                        csNongapRdRef[ ixgapalnGfend, 0 ] - 1 )

    ################

    lixgrngMatch = ixgapalnGfstart + arrayToNonZeroContigRanges( 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 0 ] == 1 ) & 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 1 ] == 1 ) )
    lixgrngMatch = np.c_[ lixgrngMatch, np.repeat(0, lixgrngMatch.shape[0]) ]

    lixgrngInsToRef = ixgapalnGfstart + arrayToNonZeroContigRanges( 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 0 ] == 1 ) & 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 1 ] == 0 ) )
    lixgrngInsToRef = np.c_[ lixgrngInsToRef, np.repeat(1, lixgrngInsToRef.shape[0]) ]

    lixgrngDelFromRef = ixgapalnGfstart + arrayToNonZeroContigRanges( 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 0 ] == 0 ) & 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 1 ] == 1 ) )
    lixgrngDelFromRef = np.c_[ lixgrngDelFromRef, np.repeat(2, lixgrngDelFromRef.shape[0]) ]

    lixgrngSoftMask = ixgapalnGfstart + arrayToNonZeroContigRanges( 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 0 ] == -1 ) & 
        ( lNongapRdRef[ ixgapalnGfstart:ixgapalnGfend+1, 1 ] == 0 ) )
    lixgrngSoftMask = np.c_[ lixgrngSoftMask, np.repeat(4, lixgrngSoftMask.shape[0]) ]

    lixgrngChunks = np.r_[ lixgrngMatch,lixgrngInsToRef,lixgrngDelFromRef,lixgrngSoftMask ]
    lixgrngChunks = lixgrngChunks[ lixgrngChunks[:,0].argsort(), : ]

    newCig = []

    # must start in a match
    assert lixgrngChunks[ 0, 2 ] == 0, 'CIGAR did not start in match'

    for ichunk in range(lixgrngChunks.shape[0]):

        newCigCode = lixgrngChunks[ichunk,2]

        # the length in the gapped alignment can be used here
        newCigLen = lixgrngChunks[ ichunk, 1 ] - lixgrngChunks[ ichunk, 0 ] + 1

        newCig += [ (newCigCode, newCigLen) ]

    # DEBUG_OLD_REF_START = rd.reference_start
    # DEBUG_OLD_REF_END = rd.reference_end-1

    rd.query_sequence = rds[ rnginclQueryNew[0]:rnginclQueryNew[1]+1 ]
    rd.reference_start = cozrefNewStart
    rd.query_qualities = rdq[ rnginclQueryNew[0]:rnginclQueryNew[1]+1 ]

    rd.cigartuples = newCig 

    # if DEBUG_OLD_REF_START <= gapFillStart and DEBUG_OLD_REF_END >= gapFillEnd:
    #     print DEBUG_OLD_REF_START,DEBUG_OLD_REF_END,gapFillStart,gapFillEnd,rd.reference_start,rd.reference_end-1

    return rd


def main():               
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--inBam',default=None,dest='inBam')
    opts.add_argument('--outBam',default=None,dest='outBam')
    
    opts.add_argument('--inMipTable', dest='inMipTable')

    opts.add_argument('--debugNonCythonized', default=True, action='store_false', dest='useCythonized')

    o = opts.parse_args()

    ######################################################
    # load MIP information.

    tblMips = pd.read_csv(o.inMipTable,sep='\t')

    tblMipsByIx = tblMips.sort_values( by=['mip_index'] ).set_index('mip_index')

    mProbeCorzGf = {}
    for i,r in tblMipsByIx.iterrows():
        mProbeCorzGf[ i ] = (r.start,r.end)

    ######################################################

    bamIn = pysam.Samfile( o.inBam, 'rb' )
    bamOut = pysam.Samfile( o.outBam, 'wb', template=bamIn )
    
    ######################################################

    ctr=0
    useCythonized = o.useCythonized

    while True:

        #debug
        # if ctr == 100000:
        #     break

        if ctr%50000==0:
            sys.stderr.write('%d..'%ctr)
            sys.stderr.flush()

        try:
            r=next(bamIn)
            ctr+=1
        except StopIteration:
            break

        tags = dict(r.tags)

        ipb = tags[TAG_PROBE_ID]
        
        concordance = tags[TAG_READ_CONCORDANT]

        if ipb >= 0 and concordance==0:

            if 'MD' in tags:
                del tags['MD']
            if 'NM' in tags:
                del tags['NM']

            # accessing single rows from pandas tables 
            # inside loops is slow...

            # probe = tblMipsByIx.loc[ipb]
            # corzGapfill = (probe.start, probe.end)
            corzGapfill = mProbeCorzGf[ ipb ]

            if useCythonized:
                # clip, cythonized 
                was_handled = clipReadToGapfill_pyx( 
                    corzGapfill[0], corzGapfill[1], r )
            else:
                was_handled = clipReadToGapfill( corzGapfill[0], corzGapfill[1], r )


            if was_handled:
                r.tags = [ (k,tags[k]) for k in tags ]
                bamOut.write(r)
            
            # debug:
            assert r.infer_query_length() == len(r.query_sequence), (str(r), corzGapfill)
            assert r.cigartuples[0][0] in (0,1), r 
            assert r.cigartuples[-1][0] in (0,1), r




                    
if __name__ == '__main__': 
    main()