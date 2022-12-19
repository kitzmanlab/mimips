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
# def split_coverage(x):
#     """
#     Split a coverage file created using bedtools coverage -hist -- which will
#     have trailing "all" hist lines -- into 1) a BedTool object with valid BED
#     lines and 2) a pandas DataFrame of all coverage, parsed from the trailing
#     "all" lines.

#     `x` can be a filename or a BedTool instance.
#     """
#     if isinstance(x, basestring):
#         fn = x
#     else:
#         fn = x.fn
 
#     f = open(fn)
#     hist_lines = []
 
#     def gen():
#         """
#         Generator that yields only valid BED lines and then stops.

#         The first "all" line is appended to hist_lines.
#         """
#         while True:
#             line = next(f)
#             toks = line.strip().split('\t')
#             if toks[0] == 'all':
#                 # Don't need to include the "all" text in the first field.
#                 hist_lines.append(toks[1:])
#                 break
#             yield pbt.create_interval_from_list(toks)
 
#     # Create a BedTool from the generator, and call `saveas` to consume the
#     # generator.  We need this so that the file pointer will stop right after
#     # the first "all" line.
#     b = pbt.BedTool(gen()).saveas()
 
#     # Pick up where we left off in the open file, appending to hist_lines.
#     while True:
#         try:
#             line = next(f)
#         except StopIteration:
#             break
#         hist_lines.append(line.strip().split('\t')[1:])
 
#     df = pd.DataFrame(
#         hist_lines,
#         columns=['depth', 'count', 'size', 'percent'])

#     return b, df

def computeCoverage( fnBed, fnBam, genomeChromSizes, libname, extendReadsBy=0, groupTargsByName=False ):

    btTargets = pbt.BedTool( fnBed )

    if extendReadsBy == 0:
        btReadsIn = pbt.BedTool( fnBam )
        btReads = btReadsIn
    else:
        # btReadsIn = pbt.BedTool( fnBam ).bam_to_bed( ).cut( [0,1,2] )
        # need to capture strand info

        # bam_to_bed chokes on empty bam files... so count first...
        btReadsIn = pbt.BedTool( fnBam )

        if len(btReadsIn)>0:
            btReadsIn = pbt.BedTool( fnBam ).bam_to_bed( ).cut( [0,1,2,3,4,5] )        
            btReads = btReadsIn.slop( g=genomeChromSizes, l=0, r=extendReadsBy, s=True )
        else:
            btReads = btReadsIn


    colsOutOverallHisto = ['depth','count','size','percent','libname','cumulative_percent']
    colsOutPerTgtSummary = ['tgtchrom','tgtstart','tgtend','mean_cvg','median_cvg','std_cvg','mad_cvg']
    colsOutPerBase = ['tgtchrom','tgtstart','tgtend','pos','cvg']

    tblOverallHisto = pd.DataFrame(
        OrderedDict( [ (c,[]) for c in colsOutOverallHisto ] ))

    # if len(btReadsIn) == 0:

    #     tblPerBase = pd.DataFrame(
    #     OrderedDict( [ (c,[]) for c in colsOutPerBase] ))

    #     tblPerTargetSummary = OrderedDict( [ (c,[]) for c in colsOutPerTgtSummary] )

    #     for iv in btTargets:
    #         tblPerTargetSummary['tgtchrom'].append(iv.chrom)
    #         tblPerTargetSummary['tgtstart'].append(iv.start)
    #         tblPerTargetSummary['tgtend'].append(iv.end)
    #         tblPerTargetSummary['mean_cvg'].append(0)
    #         tblPerTargetSummary['median_cvg'].append(0)
    #         tblPerTargetSummary['std_cvg'].append(0)
    #         tblPerTargetSummary['mad_cvg'].append(0)

    #     tblPerTargetSummary = pd.DataFrame(tblPerTargetSummary)

    #     return tblOverallHisto, tblPerTargetSummary, tblPerBase

    tblPerBase = btTargets.coverage( btReads, d=True ).to_dataframe()
    # tblPerBase = tblPerBase[ tblPerBase.columns[:5] ]

    # import pdb ; pdb.set_trace()

    if 'name' not in btTargets.to_dataframe().columns:
        tblPerBase.columns = [ 'tgtchrom','tgtstart','tgtend','pos','cvg' ]
    else:
        tblPerBase.columns = [ 'tgtchrom','tgtstart','tgtend','tgtname','pos','cvg' ]

    tblPerBase['tgtchrom']=tblPerBase['tgtchrom'].astype('str')
    for col in ('tgtstart','tgtend','pos','cvg'):
        tblPerBase[col]=tblPerBase[col].astype('int32')

    tblPerBase['pos']=tblPerBase['pos'] + tblPerBase['tgtstart']  #-> 1 based

    maxCvg = tblPerBase['cvg'].max()

    arDepths = np.arange( 0, maxCvg+1 )

    histoCvg, _ = np.histogram( tblPerBase['cvg'], arDepths )

    tblOverallHisto['depth'] = arDepths[:-1]
    tblOverallHisto['count'] = histoCvg
    tblOverallHisto['size'] = tblPerBase.shape[0]
    tblOverallHisto['percent'] = 100*(tblOverallHisto['count']/float(tblOverallHisto['count'].sum()))
    tblOverallHisto['libname'] = libname
    tblOverallHisto['cumulative_percent'] = 100 - np.r_[ 0., tblOverallHisto['percent'][:-1].cumsum() ]

    # aggregate by target
    if not groupTargsByName:
        bytgt_summary = tblPerBase.groupby( ['tgtchrom','tgtstart','tgtend'] )['cvg'].agg( [np.mean, np.median, np.std, lambda x:np.median(np.abs(x-np.median(x)))] )
    else:
        bytgt_summary = tblPerBase.groupby( ['tgtname'] )['cvg'].agg( [np.mean, np.median, np.std, lambda x:np.median(np.abs(x-np.median(x)))] )

    bytgt_summary.columns = ['mean_cvg','median_cvg','std_cvg','mad_cvg']
    bytgt_summary = bytgt_summary.reset_index()

    return tblOverallHisto, bytgt_summary, tblPerBase


    # for cvgThresh in lThresholds:


    # colsOutPerTgtHisto = ['target_chrom','target_start','target_end','depth','nbp_at_depth','target_total_bp','target_frac_bp']
    # colsOutPerBase = ['tgtchrom','tgtstart','tgtend','pos','cvg']




    # # note - this does not handle the case where the read would be extended beyond the gap fill start/end.
    # if extendReadsBy>0:
    #     btReads = btReadsIn.slop( g=genomeChromSizes, l=0, r=extendReadsBy, s=True )
    # else:
    #     btReads = btReadsIn

    # btCoverage = btReads.coverage( btTargets, d=True )

    # btPerTargetCvg, tblOverallHisto = split_coverage( btCoverage )

    # tblPerTargetHisto = \
    #     pd.DataFrame(
    #         [ list(r) for r in btPerTargetCvg ],
    #         columns=colsOutPerTgtHisto )

    # for k in ['target_start','target_end','depth','nbp_at_depth','target_total_bp']: tblPerTargetHisto.dtypes[k]='int32'
    # for k in ['target_frac_bp']: tblPerTargetHisto.dtypes[k]='float64'

    # tblPerTargetHisto = tblPerTargetHisto.sort(['target_chrom','target_start','target_end'])
    # # tblPerTargetHisto = tblPerTargetHisto.set_index()
    # perTargetSummary = tblPerTargetHisto.groupby(['target_chrom','target_start','target_end']).apply( 
    #         lambda g:np.dot( g['target_frac_bp'].astype('f') , g['depth'].astype('f')  ) )

    # tblPerTargetSummary = pd.DataFrame( {'mean_depth':perTargetSummary} )

    # lCvgThreshCols=[]

    # for cvgThresh in lThresholds:
    #     tblAtCurThresh = tblPerTargetHisto.groupby(['target_chrom','target_start','target_end']).apply( 
    #             lambda g: (g['target_frac_bp'].astype('f')[ g[ 'depth' ].astype('f')>=cvgThresh ]).sum() )

    #     colname = 'frac_gte_cvg_%d'%cvgThresh
    #     lCvgThreshCols.append(colname)

    #     tblPerTargetSummary.ix[ tblAtCurThresh.index, colname ] = tblAtCurThresh


    # # print tblPerTargetSummary.columns

    # tblPerTargetSummary = tblPerTargetSummary.reset_index()
    # tblPerTargetSummary['libname'] = libname


    # tblOverallHisto['libname']=libname
    # tblOverallHisto['percent'] = tblOverallHisto['percent'].astype('f')
    # tblOverallHisto['cumulative_percent'] = 1.0 - np.array(tblOverallHisto['percent']).cumsum()
    
    # return tblOverallHisto, tblPerTargetSummary

def main():          
    
    opts = argparse.ArgumentParser()
    
    opts.add_argument('--inBam',default=None,dest='inBam')
    
    opts.add_argument('--extendReadsBy',type=int,default=0,dest='extendReadsBy')
    opts.add_argument('--genomeChromSizes',default=None,dest='genomeChromSizes')

    opts.add_argument('--libname',default=None,dest='libname')
    opts.add_argument('--inBedTargets',default=None,dest='inBedTargets')

    opts.add_argument('--groupTargsByName',default=False,action='store_true',dest='groupTargsByName')

    opts.add_argument('--outUnifTbl', dest='outUnifTbl')
    opts.add_argument('--outPerTgtTbl', dest='outPerTgtTbl')
    opts.add_argument('--outPerBaseTbl', dest='outPerBaseTbl')

    o = opts.parse_args()

    ######################################################

    tblOverallHisto, tblPerTargetSummary, tblPerBase = \
        computeCoverage( o.inBedTargets,
                         o.inBam,
                         o.genomeChromSizes,
                         o.libname,
                         o.extendReadsBy,
                         o.groupTargsByName  )    
    
    tblPerTargetSummary.to_csv(o.outPerTgtTbl,index=False,sep='\t')
    
    tblOverallHisto.to_csv(o.outUnifTbl,index=False,sep='\t')

    if o.outPerBaseTbl is not None:
        tblPerBase.to_csv(o.outPerBaseTbl, index=False, sep='\t')


if __name__ == '__main__':
    main()