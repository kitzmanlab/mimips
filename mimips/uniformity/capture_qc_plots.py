from __future__ import division
from __future__ import print_function
#from builtins import str  # this breaks all sorts of other scripts which test type(x)==str
from builtins import range
import sys
import argparse
from collections import OrderedDict

from matplotlib.colors import LogNorm


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
            
        tblout.loc[cur_r,cur_c]=cur_v

    
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


    opts.add_argument('--autoscale_readsbywell',default=False,action='store_true',dest='autoscale_readsbywell')

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

        in_qcsummary = pd.merge( _in_qcsummary, platekey, how='inner', on='libname', suffixes=['_qc',''] )

        if in_qcsummary.shape[0] < _in_qcsummary.shape[0]:
            print( 'WARNING: %d libraries in QC summary are not found in the plate key file: %s'%( 
                _in_qcsummary.shape[0] - in_qcsummary.shape[0],
                '\n'.join( list(set(_in_qcsummary['libname']).difference( set(in_qcsummary['libname']))) )
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
    # by captureset:  % target coverage @ thresh, vs num pairs input
    #    

    pointcolors = sns.color_palette("Set2", 8)

    colLab = { 'num_pairs_input':'# raw input reads',
               'num_pairs_dedup_mapped': '# dedupd mapped reads'   }

    for capset, bycapset in in_qcsummary.groupby('captureset'):

        for cvg_col in ['num_pairs_input','num_pairs_dedup_mapped']:

            f,ax = plt.subplots(1,1)

            icol = 0
            for cvgThresh in cvg_thresholds:
                bycapset.plot.scatter( x=cvg_col, y='dedupdcvg_gte_%d'%cvgThresh, label='>=%dX'%cvgThresh, ax=ax,
                    color=pointcolors[icol%len(pointcolors)] )
                icol += 1

            plt.title( capset + ', targets >= %d X (dedups) vs %s'%( cvgThresh, colLab[cvg_col] ) )


            # Shrink current axis by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            plt.ylim(0,100)


            plt.savefig( '%s_coverage_vs_%s.pdf'%(o.out_base+'_'+capset, cvg_col) )



    ################################################
    #
    # pooling balance plot, plate-wise
    #

    if o.plot_by_plate:

        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'num_pairs_input', dtype=float )

        for plate in mPlateTbl:

            hwka={}

            if o.autoscale_readsbywell:
                hwka={'annot':True, 'fmt':'.1e' }
            else:
                hwka={'annot':True, 'fmt':'.1e', 'vmin':5e5, 'vmax':3e6 }

            f,ax = plotByWellTbl(  mPlateTbl[plate],
                                   heatmap_kwargs=hwka,
                                   annot_kwargs={'size':5.}  )

            plt.title('Capture set %s, # read pairs by position'%plate)

            plt.tight_layout()

            f.savefig('%s_%s.by_well.num_readpair_perlib.pdf'%(o.out_base, plate))


    ################################################
    #
    # pooling balance plot, plate-wise
    #

    if o.plot_by_plate:

        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'num_pairs_input', dtype=float )

        for plate in mPlateTbl:
            mPlateTbl[plate]=1e8*mPlateTbl[plate]/float(mPlateTbl[plate].sum().sum())

        for plate in mPlateTbl:

            hwka={}

            hwka={'annot':True, 'fmt':'.1e', 'vmin':1e5, 'vmax':1.5e6 , 
             'norm':LogNorm(vmin=2e5,vmax=1.5e6), 
             'cbar_kws':{
                'ticks':[1e5,2.5e5,5e5,1e6,1.5e6],
                'format':'%.2e'
                },
             'cmap':'viridis'}

            f,ax = plotByWellTbl(  mPlateTbl[plate],
                                   heatmap_kwargs=hwka,
                                   annot_kwargs={'size':5.}  )

            plt.title('Capture set %s, projected # reads pairs per 100 million'%plate)

            plt.tight_layout()

            f.savefig('%s_%s.by_well.readsperHmillion.pdf'%(o.out_base, plate))


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

                f.savefig('%s_%s.by_well.pcttarget_%dxcvg_perlib.pdf'%(o.out_base, plate, cvgThresh))


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

            f.savefig('%s_%s.by_well.frac_pairs_mapped.pdf'%(o.out_base, plate))



        mPlateTbl = makeByWellTbl( in_qcsummary, o.col_platewell, o.col_platename, 'postFiltNondupMapRate', dtype=float  )

        for plate in mPlateTbl:

            f,ax = plotByWellTbl(  100.*mPlateTbl[plate],
                                   heatmap_kwargs={'annot':True, 'fmt':'.1f', 
                                                   'cmap':sns.light_palette("red",as_cmap=True),
                                                   'vmin':30.,'vmax':90.},
                                   annot_kwargs={'size':7.}  )

            plt.title('Capture set %s, %% properly mapped paired and nonduplicate reads'%plate)

            plt.tight_layout()

            f.savefig('%s_%s.by_well.postFiltNondupMapRate.pdf'%(o.out_base, plate))


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

            f.savefig('%s_%s.by_well.duprate.pdf'%(o.out_base, plate))


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

            f.savefig('%s_%s.by_well.est_lib_size.pdf'%(o.out_base, plate))





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

