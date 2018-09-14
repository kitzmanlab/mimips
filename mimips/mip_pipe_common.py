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


def coordConversion( typeFrom='[01]', typeTo='[00]' ):
    intvCorrectionFrom=[0,0]
    intvCorrectionFrom[0] = ( typeFrom[0]=='(' and 1 or
                          typeFrom[0]=='[' and 0 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 
    intvCorrectionFrom[1] = ( typeFrom[3]==')' and 0 or
                          typeFrom[3]==']' and 1 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 

    intvCorrectionTo=[0,0]
    intvCorrectionTo[0] = ( typeTo[0]=='(' and 1 or
                          typeTo[0]=='[' and 0 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 
    intvCorrectionTo[1] = ( typeTo[3]==')' and 0 or
                          typeTo[3]==']' and 1 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 

    return ( intvCorrectionFrom[0]-intvCorrectionTo[0], 
             intvCorrectionFrom[1]-intvCorrectionTo[1] )

def data_frame_to_pbt_nameidx( df, 
                               col_chrom='chrom', 
                               col_start='start', 
                               col_end='end',
                               col_idx=None,
                               use_index=True,
                               coords='[00]' ):

    import pybedtools as pbt

    if col_idx is None and use_index:
        df_tobt = df[ [ col_chrom, col_start, col_end ] ].copy()
        df_tobt['name'] = df.index
    elif col_idx is None and not use_index:
        df_tobt = df[ [ col_chrom, col_start, col_end ] ].copy()
        df_tobt['name'] = np.arange( df.shape[0] )
    else:
        df_tobt = df[ [ col_chrom, col_start, col_end, col_idx ] ].copy()

    df_tobt.columns = ['chrom','start','end','name']
        
    coord_offset = coordConversion( coords, '[00)' )
    if coord_offset!=(0,0):
        df_tobt['start']+=coord_offset[0]
        df_tobt['end']+=coord_offset[1]

    outbt = pbt.BedTool.from_dataframe( df_tobt )

    return outbt


# TODO add optional param to this and data_frame_to_pbt_nameidx to keep strand info
def bed_intersect_dataframe( bed,
                             df_table, 
                             col_chrom='chrom', 
                             col_start='start', 
                             col_end='end',
                             col_idx=None,
                             use_index=True,
                             coords='[00]',
                             *intersect_args ):

    import pybedtools as pbt
    import numpy as np
  
    bt_tbl = data_frame_to_pbt_nameidx( df_table, col_chrom, col_start, col_end, col_idx, use_index, coords )

    bt_tbl_ovl = bt_tbl.intersect( bed, *intersect_args )

    # get indices from bed intersect results
    # cast back to proper type to index this table.
    lidcs = np.unique( np.array( [ iv.fields[3] for iv in bt_tbl_ovl ] ).astype( df_table.index.dtype ) )

    df_tbl_out = df_table.ix[ lidcs ]

    return df_tbl_out


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


