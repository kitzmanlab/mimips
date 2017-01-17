from collections import OrderedDict 
import numpy as np
import gzip

tabnl=lambda l:'\t'.join([str(x) for x in l])+'\n'

def zopen(fn,mode):
    return gzip.GzipFile(fn,mode) if fn.endswith('gz') else open(fn,mode)

# for an array, return inclusive ranges which are continuously
# equal to a single value (excluding ranges equal to zeroDefVal)    
def arrayToNonZeroContigRanges( A, zeroDefVal=0 ):
    dA=np.diff(A)
    N=A.shape[0]
    nzdA,=np.nonzero(dA)
    liRangesAll=np.c_[ np.r_[0,nzdA+1],
                       np.r_[nzdA,N-1] ]
    iiR = A.take(liRangesAll[:,0]) != zeroDefVal
    return liRangesAll[iiR,:]


#
def parseContigLengths(fnContigLengths):
    mCn=OrderedDict()
    filContigLengths = open(fnContigLengths,'r')
    l=filContigLengths.readline()
    if l[0]=='@':
        while len(l)>0:
            l=l.replace('\n','').split('\t')
            if l[0]=='@SQ':
                m={}
                for kv in l[1:]: m[kv.split(':')[0]]=kv.split(':')[1]
                if 'SN' in m and 'LN' in m:
                    mCn[m['SN']]=int(m['LN'])
            l=filContigLengths.readline()            
    else:
        while len(l)>0:
            l=l.replace('\n','').split('\t')
            mCn[l[0]]=int(l[1])
            l=filContigLengths.readline()
    return mCn

"""
	Tags added to bam file

	YT - molecular ID

	ZE - mapping to MIP based upon ext probe
            -1 means fwd read unmapped
            -2 means did not match any ext probe within tolerance
	ZF - number of ext probes within tol


	ZL - mapping to MIP based upon lig probe
            -1 means rev read unmapped
            -2 means did not match any lig probe within tolerance
	ZM - number of lig probes within tol


	ZP - mapping to MIP based upon ext and lig probe
            -1 means either read unmapped
            -2 means either did not match within tolerance
            -3 means reads on different chroms
            -4 means matched arms are not paired together in MIPS tbl
	ZQ - number of MIP probes within tol

	ZC - pair and insert size is concordant with MIP probe and capture of insert 
            0 OK
            1 read1 not cover gap fill enough
            2 read2 not cover gap fill enough
            4 read1 not on forward strand
            8 read2 not on reverse strand


    numReadsOK = YP>=0 and YC>=0
    numReadsAcap = YP>=0 and YC in (1,2)
    numReadsFilteredOther = (total)-numReadsOK-numReadsAcap
"""


TAG_EXTARM_PROBE_ID="ZE"
TAG_EXTARM_PROBE_NCANDS="ZF"

TAG_LIGARM_PROBE_ID="ZL"
TAG_LIGARM_PROBE_NCANDS="ZM"

TAG_PROBE_ID="ZP"
TAG_PROBE_NCANDS="ZQ"

TAG_READ_CONCORDANT="ZC"

mip_bam_tags = (\
	TAG_EXTARM_PROBE_ID,
	TAG_EXTARM_PROBE_NCANDS,
	TAG_LIGARM_PROBE_ID,
	TAG_LIGARM_PROBE_NCANDS,
	TAG_PROBE_ID,
	TAG_PROBE_NCANDS,
	TAG_READ_CONCORDANT)

mip_tag_types = {
	TAG_EXTARM_PROBE_ID: 'i',
	TAG_EXTARM_PROBE_NCANDS: 'i',
	TAG_LIGARM_PROBE_ID: 'i',
	TAG_LIGARM_PROBE_NCANDS: 'i',
	TAG_PROBE_ID: 'i',
	TAG_PROBE_NCANDS: 'i',
	TAG_READ_CONCORDANT: 'i' }


