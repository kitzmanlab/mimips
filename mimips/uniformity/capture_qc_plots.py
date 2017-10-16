from __future__ import division
from __future__ import print_function
#from builtins import str  # this breaks all sorts of other scripts which test type(x)==str
from builtins import range
from past.utils import old_div
import sys
import argparse
from collections import OrderedDict

import os.path

import pandas as pd

import numpy as np

import pylab as plt
import seaborn as sns

# returns dict of captureset -> perwelltbl
def makeByWellTbl( sqctbl,   
                   wellcol='well',  
                   platecol='platename',
                   valcol=None, 
                   fxnRecToValue=None,
                   dtype=None):
    
    res = {}
    
    for capset, capsettbl in sqctbl.groupby(platecol):
        res[capset] = makeWellTblSingleCapset( capsettbl, wellcol, valcol, fxnRecToValue, dtype )
    
    return res

def makeWellTblSingleCapset( sqctbl,   
                               wellcoll='well',  
                               valcol=None, 
                             fxnRecToValue=None, 
                             dtype=None ):
    
    platecols=[str(x) for x in range(1,13)]
    platerows='ABCDEFGH'
    
    # if dtype is None:
    #     tblout = OrderedDict( [ (str(c), OrderedDict( [ (r,None) for r in platerows ])) for c in platecols ] )
    # else:
    #     tblout = OrderedDict( [ (str(c), OrderedDict( [ (r,dtype()) for r in platerows ])) for c in platecols ] )

    tblout = OrderedDict( [ (str(c), OrderedDict( [ (r,None) for r in platerows ])) for c in platecols ] )

    tblout = pd.DataFrame( tblout )
    for c in tblout.columns:
        tblout[c] = tblout[c].astype(dtype)

    for _,r in sqctbl.iterrows():
        cur_rc=r[wellcoll]
        cur_r = cur_rc[0]
        cur_c = cur_rc[1:]
        
        if valcol is not None:
            cur_v = r[valcol]
        else:
            cur_v = fxnRecToValue( r )
        
        if dtype is not None:
            cur_v = dtype(cur_v)
            
        tblout.ix[cur_r,cur_c]=cur_v

    
    return tblout


def plotByWellTbl( bywelltbl,
                   title='',
                   heatmap_kwargs={},
                   annot_kwargs={} ):
        
    f,ax=plt.subplots(1,1)
    
    sns.heatmap( data=bywelltbl, ax=ax, annot_kws=annot_kwargs, **heatmap_kwargs )
    
    plt.title(title)    
    
    return f,ax
    


if __name__ == '__main__':                
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--fxn_row_to_captureset',default='lambda r:"ALL"', dest='fxn_row_to_captureset',
        help='specify a lambda function to map a row to a captureset; output plots will be grouped by capture set')

    opts.add_argument('--in_qcsummary',dest='in_qcsummary')
    opts.add_argument('--out_base',default='out',dest='out_base')
    opts.add_argument('--cvg_thresholds',default='4,8,20,100',dest='cvg_thresholds')

    # for plotting by capture plate
    # 
    #  this key must define plate name (in column given by --col_platename), plate well (--col_platewell) and "libname"

    opts.add_argument('--plot_by_plate',default=False,action='store_true',dest='plot_by_plate')

    opts.add_argument('--in_platekey', dest='in_platekey')
    opts.add_argument('--col_platename', dest='col_platename')
    opts.add_argument('--col_platewell', dest='col_platewell')

    o = opts.parse_args()

    fxn_row_to_captureset = eval(o.fxn_row_to_captureset)

    cvg_thresholds = [int(x) for x in o.cvg_thresholds.split(',')]

    in_qcsummary = pd.read_table( o.in_qcsummary )

    for i, r in in_qcsummary.iterrows():
        in_qcsummary.loc[ i, 'captureset' ] = fxn_row_to_captureset( r )

    in_qcsummary['postFiltNondupMapRate'] = (1.+in_qcsummary['num_pairs_dedup_mapped']) / (1.+in_qcsummary['num_pairs_input'])

    in_qcsummary = in_qcsummary.sort_values( by=['captureset','libname'] )

    if o.plot_by_plate:
        platekey = pd.read_table( o.in_platekey )
        assert o.col_platename in platekey.columns, "must provide --col_platename to define which column in platekey specifies the plate for each library"
        assert o.col_platewell in platekey.columns, "must provide --col_platewell to define which column in platekey specifies the plate for each library"

        platekey = platekey[ ['libname',o.col_platename,o.col_platewell] ]

        _in_qcsummary = in_qcsummary.copy()
        in_qcsummary = pd.merge( _in_qcsummary, platekey, how='inner', on='libname'  )

        if in_qcsummary.shape[0] < _in_qcsummary.shape[0]:
            print( 'WARNING: %d libraries in QC summary are not found in the plate key file: %s'%( 
                _in_qcsummary.shape[0] - in_qcsummary.shape[0],
                '\n'.join( list(set(_in_qcsummary['libname']).difference( set(in_qcsummary['libnmae']))) )
            ) )

    ################################################
    #
    # pooling balance plot
    #

    for capset, bycapset in in_qcsummary.groupby('captureset'):

        # compute MAD( #input reads ) over libraries
        liNonemptyLibs = ~(bycapset['libname'].str.lower().str.contains('emptycap') | \
                           bycapset['libname'].str.lower().str.contains('ntc') )
        
        medianNumPairs = np.median( bycapset['num_pairs_input'][liNonemptyLibs] )

        madNumPairs = np.median( np.abs( bycapset['num_pairs_input'][liNonemptyLibs] - medianNumPairs ) )
        pctmadNumPairs = 100.*madNumPairs/float(medianNumPairs)

        print('MAD(#pairs)=%0.2e, %.2f%% of median(#pairs)'%( madNumPairs,  pctmadNumPairs ))

        bycapsetSortByNrinput = bycapset.sort_values( by='num_pairs_input' ).set_index('libname')

        bycapsetSortByNrinput[ ['num_pairs_input','num_pairs_mapped','num_pairs_dedup_mapped'] ].plot.bar(stacked=False, lw=0, figsize=(12,6))

        plt.ylabel('# read pairs per library')

        plt.title('Pooling balance\nMAD(#input pairs)=%0.2e, %.2f%% of median(#input pairs)'%( madNumPairs,  pctmadNumPairs ))

        plt.tight_layout()

        plt.savefig('%s_num_readpair_perlib.pdf'%(o.out_base+'_'+capset))


    ################################################
    #
    # coverage vs num pairs input
    #

    # for cvgThresh in cvg_thresholds:

    #     f,ax = plt.subplots(1)

    #     plt.scatter( )


    #     plt.title('Capture set %s, log10(estimated library size)'%plate)

    #     plt.tight_layout()

    #     f.savefig('%s.%s.by_well.est_lib_size.pdf'%(o.out_base, plate))




    ################################################
    #
    # pooling balance plot, plate-wise
    #

    if o.plot_by_plate:

        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'num_pairs_input', dtype=float )

        for plate in mPlateTbl:

            f,ax = plotByWellTbl(  mPlateTbl[plate],
                                   heatmap_kwargs={'annot':True, 'fmt':'.1e',
                                   'vmin':5e5, 'vmax':3e6 },
                                   annot_kwargs={'size':7.}  )

            plt.title('Capture set %s, # read pairs by position'%plate)

            plt.tight_layout()

            f.savefig('%s%s.by_well.num_readpair_perlib.pdf'%(o.out_base, plate))


    ################################################
    #
    # dedupdcvg_gte_{x} platewise
    #


    if o.plot_by_plate:

        for cvgThresh in cvg_thresholds:

            mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'dedupdcvg_gte_%d'%cvgThresh, dtype=float  )

            for plate in mPlateTbl:

                f,ax = plotByWellTbl(  mPlateTbl[plate],
                                       heatmap_kwargs={'annot':True, 'fmt':'.1f', 'cmap':'coolwarm','vmin':80,'vmax':100},
                                       annot_kwargs={'size':7.}  )

                plt.title('Capture set %s, %% target bases >= %dX'%(plate,cvgThresh))

                plt.tight_layout()

                f.savefig('%s%s.by_well.pcttarget_%dxcvg_perlib.pdf'%(o.out_base, plate, cvgThresh))


    ################################################
    #
    # map rate platewise
    #


    if o.plot_by_plate:

        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'frac_pairs_mapped', dtype=float  )

        for plate in mPlateTbl:

            f,ax = plotByWellTbl(  100.*mPlateTbl[plate],
                                   heatmap_kwargs={'annot':True, 'fmt':'.1f', 
                                                   'cmap':sns.light_palette("red",as_cmap=True),
                                                   'vmin':30.,'vmax':90.},
                                   annot_kwargs={'size':7.}  )

            plt.title('Capture set %s, %% properly mapped paired reads'%plate)

            plt.tight_layout()

            f.savefig('%s%s.by_well.frac_pairs_mapped.pdf'%(o.out_base, plate))



        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'postFiltNondupMapRate', dtype=float  )

        for plate in mPlateTbl:

            f,ax = plotByWellTbl(  100.*mPlateTbl[plate],
                                   heatmap_kwargs={'annot':True, 'fmt':'.1f', 
                                                   'cmap':sns.light_palette("red",as_cmap=True),
                                                   'vmin':30.,'vmax':90.},
                                   annot_kwargs={'size':7.}  )

            plt.title('Capture set %s, %% properly mapped paired and nonduplicate reads'%plate)

            plt.tight_layout()

            f.savefig('%s%s.by_well.postFiltNondupMapRate.pdf'%(o.out_base, plate))


    ################################################
    #
    # duprate platewise
    #


    if o.plot_by_plate:

        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'duprate', dtype=float  )

        for plate in mPlateTbl:

            f,ax = plotByWellTbl(  100.*mPlateTbl[plate],
                                   heatmap_kwargs={'annot':True, 'fmt':'.1f', 
                                                   'cmap':sns.light_palette("orange",as_cmap=True),
                                                   'vmin':0,'vmax':30},
                                   annot_kwargs={'size':7.}  )

            plt.title('Capture set %s, %% read duplication rate'%plate)

            plt.tight_layout()

            f.savefig('%s%s.by_well.duprate.pdf'%(o.out_base, plate))


    ################################################
    #
    # log(libsize) platewise
    #


    if o.plot_by_plate:


        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'est_lib_size', dtype=float  )

        for plate in mPlateTbl:

            mPlateTbl[plate] = np.log10( mPlateTbl[plate].clip(1,1e99) )

            f,ax = plotByWellTbl(  mPlateTbl[plate],
                                   heatmap_kwargs={'annot':True, 'fmt':'.1f', 
                                                   'cmap':sns.light_palette("blue",as_cmap=True),
                                                   'vmin':5,'vmax':10},
                                   annot_kwargs={'size':7.}  )

            plt.title('Capture set %s, log10(estimated library size)'%plate)

            plt.tight_layout()

            f.savefig('%s%s.by_well.est_lib_size.pdf'%(o.out_base, plate))





    #######
    # 
    # prioritize libraries for further sequencing
    #  - P(dup) ~ 1 - (cur_n_mapped_reads/2) / estd_lib_szie
    #  - 


    # todo log(libsize)


    ################################################
    #
    # cvg 8x vs # map
    #





    # could put it all in a notebook...
    # from nbformat import v4 as nbf

    # cells = [ nbf.new_code_cell( ["1+2"] ), 
    #           nbf.new_code_cell( ["2+3"] ) ]

    # nb = nbf.new_notebook( cells=cells )

    # with open('my_notebook.ipynb', 'w') as f:
    #     f.write(nbf.writes(nb))

