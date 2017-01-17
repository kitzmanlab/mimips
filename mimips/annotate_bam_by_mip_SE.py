from __future__ import division
from builtins import range
from past.utils import old_div
import sys
import argparse
from collections import defaultdict,OrderedDict

import os.path

import pysam

import pandas as pd

from mip_pipe_common import *

def matchAndTagChunkOfReads(
        chunkAlns,
        mChromLidxAlns,
        mChromPbsBoundsIdsSortByExt,
        mChromPbsBoundsIdsSortByLig,
        tblPbsByIx,
        minBpIntoGapFill,
        outerBpTolerance=2 ):
    
    Ninchunk = len(chunkAlns)
    assert len(chunkAlns)==Ninchunk

    mTagOut={}
    for tag in mip_bam_tags:
        mTagOut[tag] = -999999*np.ones( (Ninchunk,), dtype=np.int32 )


    lliCandProbes_1=[ [] for _ in range(Ninchunk) ]
    lliCandProbes_2=[ [] for _ in range(Ninchunk) ]
    
    # loop through read 1 alignments by chromosome
    for chromid in mChromLidxAlns:

        liAlns = mChromLidxAlns[chromid]

        # is read mapped?
        if chromid==-1:
            # mark all nonmapped
            mTagOut[ TAG_EXTARM_PROBE_ID ][ liAlns ] = -1
            mTagOut[ TAG_EXTARM_PROBE_NCANDS ][ liAlns ] = 0
            mTagOut[ TAG_LIGARM_PROBE_ID ][ liAlns ] = -1
            mTagOut[ TAG_LIGARM_PROBE_NCANDS ][ liAlns ] = 0

        # are there no probes on this chrom?
        elif chromid not in mChromPbsBoundsIdsSortByExt:
            # record that there were no candidates
            mTagOut[ TAG_EXTARM_PROBE_NCANDS ][ liAlns ] = 0            
            mTagOut[ TAG_LIGARM_PROBE_NCANDS ][ liAlns ] = 0            

        else:

            arPbsBoundsIdsSortByExt = mChromPbsBoundsIdsSortByExt[chromid]

            # reference_start is 0-based, inclusive 
            # so is extension probe start coord
            # -> no adjustment needed

            liipbCur = arPbsBoundsIdsSortByExt[:,0].searchsorted( 
                    [ chunkAlns[ial].reference_start - outerBpTolerance 
                      for ial in liAlns ] )

            for iialn in range(len(liAlns)):       

                iipbCur = liipbCur[iialn]

                ialn = liAlns[iialn]

                while iipbCur < arPbsBoundsIdsSortByExt.shape[0] and \
                        abs( chunkAlns[ialn].reference_start - 
                            arPbsBoundsIdsSortByExt[iipbCur,0] ) < outerBpTolerance:

                    # get probe id
                    ipbCur = arPbsBoundsIdsSortByExt[iipbCur,4]
                    lliCandProbes_1[ ialn ].append( ipbCur )
                    iipbCur+=1

                mTagOut[ TAG_EXTARM_PROBE_NCANDS ][ ialn ] = len(lliCandProbes_1[ ialn ])

            arPbsBoundsIdsSortByLig = mChromPbsBoundsIdsSortByLig[chromid]

            # reference_end is 0-based, exclusive
            # ligation arm end coord is 0-based, inclusive
            # --> use ref_end-1

            liipbCur = np.clip(
                    arPbsBoundsIdsSortByLig[:,3].searchsorted( 
                    [ chunkAlns[ial].reference_end - 1 + outerBpTolerance 
                      for ial in liAlns ], side='right' ) - 1,
                    0, 
                    arPbsBoundsIdsSortByLig.shape[0]-1 ) 

            for iialn in range(len(liAlns)):               

                iipbCur = liipbCur[iialn]

                ialn = liAlns[iialn]

                while iipbCur >= 0 and \
                      abs( chunkAlns[ialn].reference_end - 1 - 
                            arPbsBoundsIdsSortByLig[iipbCur,3] ) < outerBpTolerance :

                    # get probe id
                    ipbCur = arPbsBoundsIdsSortByLig[ iipbCur, 4 ]
                    lliCandProbes_2[ ialn ].append( ipbCur )
                    iipbCur-=1

                mTagOut[ TAG_LIGARM_PROBE_NCANDS ][ ialn ] = len(lliCandProbes_2[ ialn ])

    # loop through probes
    for ialn in range(Ninchunk):
        
        # was unmapped?
        if mTagOut[ TAG_EXTARM_PROBE_ID ][ ialn ] == -1 or\
           mTagOut[ TAG_LIGARM_PROBE_ID ][ ialn ] == -1 :

            mTagOut[TAG_PROBE_ID][ialn] = -1
            mTagOut[TAG_PROBE_NCANDS][ialn] = 0

        elif mTagOut[ TAG_EXTARM_PROBE_NCANDS ][ ialn ] == 0 or\
            mTagOut[ TAG_LIGARM_PROBE_NCANDS ][ ialn ] == 0 :

            mTagOut[TAG_PROBE_ID][ialn] = -2
            mTagOut[TAG_PROBE_NCANDS][ialn] = 0

        else:
            # mapped && can we identify a probe from both ends?
            sPbidsAgree = set( lliCandProbes_1[ialn] ).intersection( set(lliCandProbes_2[ialn]) )

            if len(sPbidsAgree) == 0:
                # no probes in agreement beteween ext and lig.

                mTagOut[TAG_PROBE_ID][ialn] = -4
                mTagOut[TAG_PROBE_NCANDS][ialn] = 0

                if len(lliCandProbes_1[ialn]) == 0:
                    mTagOut[TAG_EXTARM_PROBE_ID][ialn] = -1
                else:
                    mTagOut[TAG_EXTARM_PROBE_ID][ialn] = lliCandProbes_1[ialn][0]

                if len(lliCandProbes_2[ialn]) == 0:
                    mTagOut[TAG_LIGARM_PROBE_ID][ialn] = -1
                else:
                    mTagOut[TAG_LIGARM_PROBE_ID][ialn] = lliCandProbes_2[ialn][0]

            else:
                pbid = tuple(sPbidsAgree)[0]
                mTagOut[TAG_PROBE_ID][ialn] = pbid
                mTagOut[TAG_EXTARM_PROBE_ID][ialn] = pbid
                mTagOut[TAG_LIGARM_PROBE_ID][ialn] = pbid

                mTagOut[TAG_PROBE_NCANDS][ialn] = len(sPbidsAgree)

                # fudgeExt=tblPbsByIx.loc[ pbid ].ext_probe_start-chunkAlns_1[ialn].reference_start
                # fudgeLig=tblPbsByIx.loc[ pbid ].lig_probe_end-(chunkAlns_2[ialn].reference_end-1)
                # print fudgeExt,fudgeLig
                

    # reads are concordant when they:
    # (1) can be uniquely assigned to a probe and
    # (2) overlap the gapfill region by >= minBpIntoGapFill

    liAlnInvalidPbs=np.where(mTagOut[TAG_PROBE_NCANDS]!=1)[0]

    if liAlnInvalidPbs.shape[0]>0:
        mTagOut[TAG_READ_CONCORDANT][liAlnInvalidPbs] = 16

    liAlnValidPbs=np.where(mTagOut[TAG_PROBE_NCANDS]==1)[0]

    if liAlnValidPbs.shape[0]>0:

        arByAlnPbStartEnd = np.c_[
            tblPbsByIx.loc[ mTagOut[TAG_PROBE_ID][liAlnValidPbs] ]['start'].values,
            tblPbsByIx.loc[ mTagOut[TAG_PROBE_ID][liAlnValidPbs] ]['end'].values ]

        for iialn in range(len(liAlnValidPbs)):

            ialn = liAlnValidPbs[iialn]

            res = 0

            rd=chunkAlns[ialn]

            cozGfStart=arByAlnPbStartEnd[iialn, 0]
            cozGfEnd=arByAlnPbStartEnd[iialn, 1]

            if rd.get_overlap( cozGfStart,cozGfEnd+1 ) <= minBpIntoGapFill:
                res = res | 3

            if rd.is_reverse:
                res = res | 12

            mTagOut[TAG_READ_CONCORDANT][ialn] = res


    # assign tags
    for ialn in range(Ninchunk):

        rd=chunkAlns[ialn]
        
        for tag in mip_bam_tags:
            rd.set_tag( tag, mTagOut[tag][ialn], mip_tag_types[tag] )

    numReadsOK = (mTagOut[ TAG_READ_CONCORDANT ] == 0).sum()
    numReadsAcap = \
        ( ( mTagOut[ TAG_READ_CONCORDANT ] & 0x3 ) != 0 ).sum()
    numReadsFilteredOther =  ( ((mTagOut[ TAG_READ_CONCORDANT ] & 0x3 ) == 0) & \
                       ( mTagOut[ TAG_READ_CONCORDANT ] != 0 ) ).sum()

    return numReadsOK,numReadsAcap,numReadsFilteredOther

def main():
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--inBam',default=None,dest='inBam')
    opts.add_argument('--outBam',default=None,dest='outBam')

    opts.add_argument('--outSummary',default=None,dest='outSummary')
    
    opts.add_argument('--minBpIntoGapFill',default=20,type=int,dest='minBpIntoGapFill')

    opts.add_argument('--rg',default=None,dest='rg')
    opts.add_argument('--lib',default=None,dest='lib')  # keyline->rg
    opts.add_argument('--samp',default=None,dest='samp')  # keyline->lib
    opts.add_argument('--platform',default='ILLUMINA',dest='platform')

    opts.add_argument('--midList',default=None,dest='midList')

    opts.add_argument('--inMipTable', dest='inMipTable')

    opts.add_argument('--prefixChar', default='Z', dest='prefixChar')

    o = opts.parse_args()

    ######################################################

    mRG = [ {'ID':o.rg, 'LB':o.lib, 'SM':o.samp, 'PL':o.platform} ]

    bamIn = pysam.Samfile( o.inBam, 'rb' )

    hdrOut = bamIn.header
    hdrOut['RG']=mRG
    
    bamOut = pysam.Samfile( o.outBam, 'wb', header=hdrOut )
    
    theRg = o.rg

    filInMidList = zopen(o.midList,'r')

    numReadsOK,numReadsAcap,numReadsFilteredOther=0,0,0

    minBpIntoGapFill = o.minBpIntoGapFill

    ######################################################
    # load MIP information.

    tblMips = pd.read_csv(o.inMipTable,sep='\t')

    tblMips['chrom']=tblMips['chrom'].astype('str')

    # make by-chromosome sorted lists of positions and IDs
    mChromPbsBoundsIdsSortByExt={}
    mChromPbsBoundsIdsSortByLig={}

    # note that the coordinates loaded from the MIP design table are 0-based inclusive [00]
    for chrom, subtbl in tblMips.sort( ['chrom','ext_probe_start'] ).groupby('chrom'):

        # convert chrom -> id to match bam flines
        chromid = bamIn.gettid( chrom )

        # extension probes - track by extension probe start coord
        mChromPbsBoundsIdsSortByExt[ chromid ] = \
            np.c_[ subtbl['ext_probe_start'].values,
                   subtbl['ext_probe_end'].values,
                   subtbl['lig_probe_start'].values,
                   subtbl['lig_probe_end'].values,
                   subtbl['mip_index'].values ]
        mChromPbsBoundsIdsSortByLig[ chromid ] = mChromPbsBoundsIdsSortByExt[ chromid ].copy()

        mChromPbsBoundsIdsSortByExt[ chromid ] = \
            mChromPbsBoundsIdsSortByExt[ chromid ][
                np.argsort(mChromPbsBoundsIdsSortByExt[ chromid ][:,0]) ]

        # extension probes - track by ligation probe end coord
        mChromPbsBoundsIdsSortByLig[ chromid ] = \
            mChromPbsBoundsIdsSortByLig[ chromid ][
                np.argsort(mChromPbsBoundsIdsSortByLig[ chromid ][:,2]) ]


    tblMipsByIx = tblMips.copy().sort( ['mip_index'] ).set_index('mip_index')

    ######################################################

    chunkSize=100000

    ctr=0

    while True:

        # current chunk list of MIDs and (read1,read2)'s
        lMids,lAlns=[],[]
        
        # current chunk # reads stored
        ctrChunk=0

        # chrom --> [ list of idcs into chunk ], separately for read1 and read2
        mChromLidxCurChunkAln=defaultdict(list)
        
        # read in chunkSize lines
        while ctrChunk < chunkSize:

            ctr+=1

            if ctr%10000==0:
                sys.stderr.write('%d..'%ctr)
                sys.stderr.flush()

            # if ctr%100000==0:
            #     sys.exit(0)

            try:
                rd=next(bamIn)
            except StopIteration:
                break

            midrn=filInMidList.readline().rstrip().split('#')[0].split('/')[0].split(' ')[0][1:]
            mid = filInMidList.readline().rstrip()
            filInMidList.readline()
            filInMidList.readline()

            assert midrn == rd.qname.split('#')[0].split('/')[0].split(' ')[0]

            rd.set_tag( 'YT', mid )

            lMids.append(mid)
            lAlns.append(rd)

            chromid = rd.reference_id if not rd.is_unmapped else -1

            mChromLidxCurChunkAln[ chromid ].append(ctrChunk)

            ctrChunk+=1 

        sys.stderr.write('\nprocessing %d rows..'%ctrChunk)
        sys.stderr.flush()

        if ctrChunk==0:
            break

        curReadsOK,curReadsAcap,curReadsFilteredOther = \
            matchAndTagChunkOfReads(
                lAlns,
                mChromLidxCurChunkAln,
                mChromPbsBoundsIdsSortByExt,
                mChromPbsBoundsIdsSortByLig,
                tblMipsByIx,
                minBpIntoGapFill,
                outerBpTolerance=2 )

        numReadsOK+=curReadsOK
        numReadsAcap+=curReadsAcap
        numReadsFilteredOther+=curReadsFilteredOther

        sys.stderr.write(' running tally: %.2f%% OK, %.2f%% acap, %.2f%% other\n'%(
            100*numReadsOK/float(1+numReadsOK+numReadsAcap+numReadsFilteredOther),
            100*numReadsAcap/float(1+numReadsOK+numReadsAcap+numReadsFilteredOther),
            100*numReadsFilteredOther/float(1+numReadsOK+numReadsAcap+numReadsFilteredOther),
                ))
        sys.stderr.flush()

        for rd in lAlns:
            rd.set_tag('RG',theRg)

        for ialn in range(len(lAlns)):
            bamOut.write(lAlns[ialn])


    if o.outSummary is not None:
        tblSummary=open(o.outSummary,'w')
        tblSummary.write( tabnl( ['libname','num_reads_ok','num_reads_filter_acapture','num_reads_filter_other',
                                            'frac_reads_ok','frac_reads_filter_acapture','frac_reads_filter_other' ] ))

        fracReadsOK = old_div(numReadsOK, float( numReadsOK+numReadsFilteredOther+numReadsAcap )) if numReadsOK+numReadsFilteredOther+numReadsAcap>0 else 0.
        fracReadsAcap = old_div(numReadsAcap, float( numReadsOK+numReadsFilteredOther+numReadsAcap )) if numReadsOK+numReadsFilteredOther+numReadsAcap>0 else 0.
        fracReadsFilteredOther = old_div(numReadsFilteredOther, float( numReadsOK+numReadsFilteredOther+numReadsAcap )) if numReadsOK+numReadsFilteredOther+numReadsAcap>0 else 0.

        tblSummary.write( tabnl( [o.lib, numReadsOK, numReadsAcap, numReadsFilteredOther,
                                        fracReadsOK, fracReadsAcap, fracReadsFilteredOther ] ))

if __name__ == '__main__':                
    main()