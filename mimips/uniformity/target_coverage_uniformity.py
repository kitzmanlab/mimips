from past.builtins import basestring
import sys
import argparse
from collections import defaultdict,OrderedDict

import os.path

import pysam

import pandas as pd

import pybedtools as pbt

from ..mip_pipe_common import *


# from daler
# https://gist.github.com/daler/99cdf714d438e161dad6
def split_coverage(x):
    """
    Split a coverage file created using bedtools coverage -hist -- which will
    have trailing "all" hist lines -- into 1) a BedTool object with valid BED
    lines and 2) a pandas DataFrame of all coverage, parsed from the trailing
    "all" lines.

    `x` can be a filename or a BedTool instance.
    """
    if isinstance(x, basestring):
        fn = x
    else:
        fn = x.fn
 
    f = open(fn)
    hist_lines = []
 
    def gen():
        """
        Generator that yields only valid BED lines and then stops.

        The first "all" line is appended to hist_lines.
        """
        while True:
            line = next(f)
            toks = line.strip().split('\t')
            if toks[0] == 'all':
                # Don't need to include the "all" text in the first field.
                hist_lines.append(toks[1:])
                break
            yield pbt.create_interval_from_list(toks)
 
    # Create a BedTool from the generator, and call `saveas` to consume the
    # generator.  We need this so that the file pointer will stop right after
    # the first "all" line.
    b = pbt.BedTool(gen()).saveas()
 
    # Pick up where we left off in the open file, appending to hist_lines.
    while True:
        try:
            line = next(f)
        except StopIteration:
            break
        hist_lines.append(line.strip().split('\t')[1:])
 
    df = pd.DataFrame(
        hist_lines,
        columns=['depth', 'count', 'size', 'percent'])

    print '--------inside split_coverage-----'
    print df.shape
    print df.columns
    print df.head()
    print '-----------------------------------'


    return b, df

def computeCoverage( fnBed, fnBam, lThresholds, genomeChromSizes, libname, extendReadsBy=0 ):

    btTargets = pbt.BedTool( fnBed )

    btReadsIn = pbt.BedTool( fnBam ).bam_to_bed( ).cut( [0,1,2,3] )

    colsOutOverallHisto = ['depth','count','size','percent','libname','cumulative_percent']
    colsOutPerTgtHisto = ['target_chrom','target_start','target_end','depth','nbp_at_depth','target_total_bp','target_frac_bp']

    if len(btReadsIn) == 0:
        tblOverallHisto = pd.DataFrame(
            OrderedDict( [ (c,[]) for c in colsOutOverallHisto ] ))

        tblPerTargetHisto = pd.DataFrame(
            OrderedDict( [ (c,[]) for c in colsOutPerTgtHisto] ))

        return tblOverallHisto, tblPerTargetHisto


    # note - this does not handle the case where the read would be extended beyond the gap fill start/end.
    if extendReadsBy>0:
        btReads = btReadsIn.slop( g=genomeChromSizes, l=0, r=extendReadsBy, s=True )
    else:
        btReads = btReadsIn

    btCoverage = btReads.coverage( btTargets, hist=True )

    btPerTargetCvg, tblOverallHisto = split_coverage( btCoverage )

    print '--------- btPerTargetCvg ----------'
    print list(btPerTargetCvg[0])
    print list(colsOutPerTgtHisto)
    print '-----------------------------------'


    tblPerTargetHisto = \
        pd.DataFrame(
            [ list(r) for r in btPerTargetCvg ],
            columns=colsOutPerTgtHisto )

    for k in ['target_start','target_end','depth','nbp_at_depth','target_total_bp']: tblPerTargetHisto.dtypes[k]='int32'
    for k in ['target_frac_bp']: tblPerTargetHisto.dtypes[k]='float64'

    tblPerTargetHisto = tblPerTargetHisto.sort(['target_chrom','target_start','target_end'])
    # tblPerTargetHisto = tblPerTargetHisto.set_index()
    perTargetSummary = tblPerTargetHisto.groupby(['target_chrom','target_start','target_end']).apply( 
            lambda g:np.dot( g['target_frac_bp'].astype('f') , g['depth'].astype('f')  ) )

    tblPerTargetSummary = pd.DataFrame( {'mean_depth':perTargetSummary} )

    lCvgThreshCols=[]

    for cvgThresh in lThresholds:
        tblAtCurThresh = tblPerTargetHisto.groupby(['target_chrom','target_start','target_end']).apply( 
                lambda g: (g['target_frac_bp'].astype('f')[ g[ 'depth' ].astype('f')>=cvgThresh ]).sum() )

        colname = 'frac_gte_cvg_%d'%cvgThresh
        lCvgThreshCols.append(colname)

        tblPerTargetSummary.ix[ tblAtCurThresh.index, colname ] = tblAtCurThresh


    # print tblPerTargetSummary.columns

    tblPerTargetSummary = tblPerTargetSummary.reset_index()
    tblPerTargetSummary['libname'] = libname


    tblOverallHisto['libname']=libname
    tblOverallHisto['percent'] = tblOverallHisto['percent'].astype('f')
    tblOverallHisto['cumulative_percent'] = 1.0 - np.array(tblOverallHisto['percent']).cumsum()
    
    return tblOverallHisto, tblPerTargetSummary

def main():          
    
    opts = argparse.ArgumentParser()
    
    opts.add_argument('--inBam',default=None,dest='inBam')
    
    opts.add_argument('--extendReadsBy',type=int,default=0,dest='extendReadsBy')
    opts.add_argument('--genomeChromSizes',default=None,dest='genomeChromSizes')

    opts.add_argument('--libname',default=None,dest='libname')
    opts.add_argument('--inBedTargets',default=None,dest='inBedTargets')

    opts.add_argument('--outUnifTbl', dest='outUnifTbl')
    opts.add_argument('--outPerTgtTbl', dest='outPerTgtTbl')

    opts.add_argument('--bpThresholds', default=None, dest='bpThresholds')    

    o = opts.parse_args()

    ######################################################

    lThresholds = [] if o.bpThresholds is None else [ int(x) for x in o.bpThresholds.split(',') ]


    tblOverallHisto, tblPerTargetSummary = \
        computeCoverage( o.inBedTargets,
                         o.inBam,
                         lThresholds,
                         o.genomeChromSizes,
                         o.libname,
                         o.extendReadsBy  )    
    
    tblPerTargetSummary.to_csv(o.outPerTgtTbl,index=False,sep='\t')
    
    tblOverallHisto.to_csv(o.outUnifTbl,index=False,sep='\t')


if __name__ == '__main__':
    main()