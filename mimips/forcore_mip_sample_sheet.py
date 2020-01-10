from __future__ import print_function
from builtins import chr
from builtins import range

import Bio.Seq

import sys

import re
import os
import os.path
import argparse
from collections import defaultdict, OrderedDict

import shutil

import pandas as pd

import openpyxl

def ofsFrom(start,down=0,right=0):
    cs=start[0]
    rs=int(start[1:])
    cout=chr( ord(cs)+right )
    rout=rs+down
    return '%s%d'%( cout,rout )


def main():
    
    opts = argparse.ArgumentParser()
   
    opts.add_argument('--excelSheetIn', dest='excelSheetIn')
    
    opts.add_argument('--barcodeFasta',  default=None, dest='barcodeFasta')

    opts.add_argument('--fxnTransformBarcode_i7',  default="lambda barcin:barcin", dest='fxnTransformBarcode_i7') 
    opts.add_argument('--fxnTransformBarcode_i5',  default="lambda barcin:barcin", dest='fxnTransformBarcode_i5') 

    opts.add_argument('--coreTemplateIn', dest='coreTemplateIn')
    opts.add_argument('--coreSheetOut',  default=None, dest='coreSheetOut')
    opts.add_argument('--coreBuffer',  default=None, dest='coreBuffer')
    opts.add_argument('--coreConc',  default=None, dest='coreConc')
    opts.add_argument('--coreVol',  default=None, dest='coreVol')
    opts.add_argument('--coreFraglen',  default=None, dest='coreFraglen')
    opts.add_argument('--coreShortcode',  default=None, dest='coreShortcode')
    opts.add_argument('--coreSpecies',  default=None, dest='coreSpecies')
    opts.add_argument('--coreNotes',  default=None, dest='coreNotes')

    opts.add_argument('--libKeyOut',dest='libKeyOut')

    opts.add_argument('--extraCols',default=None,dest='extraCols')

    opts.add_argument('--fxnLibName',default='lambda r:"%s_%s"%(r.source_plate,r.source_well)',dest='fxnLibName')

    o = opts.parse_args()

    fxnLibName = eval(o.fxnLibName)

    fxnTransformBarcode_i5 = eval(o.fxnTransformBarcode_i5)
    fxnTransformBarcode_i7 = eval(o.fxnTransformBarcode_i7)

    # load the barcode sequences
    filInBarcs = open(o.barcodeFasta,'r')
    l = filInBarcs.readline()
    # if l[0]=='>':
    assert l[0]=='>'
    mBcNameSeq_i5i7={}
    while len(l)>0:
        bcname= l[1:].rstrip()
        bcseq=filInBarcs.readline().rstrip()
        l=filInBarcs.readline()
        assert bcname not in mBcNameSeq_i5i7, '%s present > once'%bcname
        mBcNameSeq_i5i7[bcname] = ( fxnTransformBarcode_i5( bcseq ), fxnTransformBarcode_i7( bcseq ) )

    wbin = openpyxl.load_workbook(filename=o.excelSheetIn)

    sheet_src_plate = wbin.get_sheet_by_name('SOURCE PLATE')    
    sheet_src_well = wbin.get_sheet_by_name('SOURCE WELL')    
    sheet_src_barcp5 = wbin.get_sheet_by_name('P5 BARCODE')    
    sheet_src_barcp7 = wbin.get_sheet_by_name('P7 BARCODE')

    # source plate:
    # make sure 1...12 from A6 to right
    for i in range(1,12):
        sheetloc=ofsFrom('A6',right=i)
        obs=str(int(sheet_src_plate[ sheetloc ].value))
        exp=str(i)
        assert obs==exp, 'SOURCE PLATE %s : expected %s but got %s'%( sheetloc, exp, obs )
    # make sure A..H giong from A7 down
    for i in range(0,8):
        sheetloc=ofsFrom('A7',down=i)
        obs=str(sheet_src_plate[ sheetloc ].value)
        exp=str( chr(ord('A')+i) )
        assert obs==exp, 'SOURCE PLATE %s : expected %s but got %s'%( sheetloc, exp, obs )

    # well plate:
    # make sure 1...12 from A6 to right
    for i in range(1,12):
        sheetloc=ofsFrom('A6',right=i)
        obs=str(int(sheet_src_well[ sheetloc ].value))
        exp=str(i)
        assert obs==exp, 'SOURCE WELL %s : expected %s but got %s'%( sheetloc, exp, obs )
    # make sure A..H giong from A7 down
    for i in range(0,8):
        sheetloc=ofsFrom('A7',down=i)
        obs=str(sheet_src_well[ sheetloc ].value)
        exp=str( chr(ord('A')+i) )
        assert obs==exp, 'SOURCE WELL %s : expected %s but got %s'%( sheetloc, exp, obs )

    # p7 barc:
    # make sure 1...12 from A7 to right
    for i in range(1,12):
        sheetloc=ofsFrom('A7',right=i)
        obs=str(int(sheet_src_barcp7[ sheetloc ].value))
        exp=str(i)
        assert obs==exp, 'P7 PRIMER %s : expected %s but got %s'%( sheetloc, exp, obs )
    # make sure A..H giong from A8 down
    for i in range(0,8):
        sheetloc=ofsFrom('A8',down=i)
        obs=str(sheet_src_barcp7[ sheetloc ].value)
        exp=str( chr(ord('A')+i) )
        assert obs==exp, 'P7 PRIMER %s : expected %s but got %s'%( sheetloc, exp, obs )

    # p5 barc:
    # make sure 1...12 from A7 to right
    for i in range(1,12):
        sheetloc=ofsFrom('A7',right=i)
        obs=str(int(sheet_src_barcp5[ sheetloc ].value))
        exp=str(i)
        assert obs==exp, 'p5 PRIMER %s : expected %s but got %s'%( sheetloc, exp, obs )
    # make sure A..H giong from A8 down
    for i in range(0,8):
        sheetloc=ofsFrom('A8',down=i)
        obs=str(sheet_src_barcp5[ sheetloc ].value)
        exp=str( chr(ord('A')+i) )
        assert obs==exp, 'p5 PRIMER %s : expected %s but got %s'%( sheetloc, exp, obs )


    # gather into DF
    df = OrderedDict()
    for col in ['well','source_plate','source_well','p7_barc_and_well','p5_barc_and_well','p7_barc','p5_barc','p7_barc_seq','p5_barc_seq']:
        df[col]=[]
    
    for j in range(0,12):
        for i in range(0,8):        
            sheetloc1=ofsFrom('B7',down=i,right=j)
            sheetloc2=ofsFrom('B8',down=i,right=j)
            # fun fact- excel does rows/cols opposite of PCR plates
            well=ofsFrom('A1',down=j,right=i)
            srcplate=sheet_src_plate[sheetloc1].value
            srcwell=sheet_src_well[sheetloc1].value
            curp5=sheet_src_barcp5[sheetloc2].value
            curp7=sheet_src_barcp7[sheetloc2].value

            srcplate='' if srcplate is None else str(srcplate).strip()
            srcwell='' if srcwell is None else str(srcwell).strip()
            curp5='' if curp5 is None else str(curp5).strip()
            curp7='' if curp7 is None else str(curp7).strip()

            # srcplate=srcplate.replace('-','').replace('_','')
            srcplate=re.subn( '[\'"$:\W@\n]', '', srcplate )[0]
            srcwell=re.subn( '[\'"$:\W@\n]', '', srcwell )[0]


            if any( [len(srcplate)==0, len(srcwell)==0, len(curp5)==0, len(curp7)==0] ):
                if not all ( [len(srcplate)==0, len(srcwell)==0, len(curp5)==0, len(curp7)==0] ):
                    print('WARNING: well %s is not empty in all sheets: %s %s %s %s'%( well, srcplate, srcwell, curp5, curp7 ))
            else:

                assert ':' in curp7, 'ERROR well %s p7 invalid barcode %s'%( well, curp7 )
                assert ':' in curp5, 'ERROR well %s p5 invalid barcode %s'%( well, curp5 )

                df['well'].append( well )
                df['source_plate'].append( srcplate )
                df['source_well'].append( srcwell )
                df['p7_barc_and_well'].append( curp7 )
                df['p5_barc_and_well'].append( curp5 )
                df['p7_barc'].append( curp7.split(':')[1] )
                df['p5_barc'].append( curp5.split(':')[1] )
                df['p7_barc_seq'].append( mBcNameSeq_i5i7[ curp7.split(':')[1] ][1] )
                df['p5_barc_seq'].append( mBcNameSeq_i5i7[ curp5.split(':')[1] ][0] )

    # gather extra cols if any
    mKvExtra={}
    if o.extraCols is not None:
        for kv in o.extraCols.split(','):
            mKvExtra[kv.split(':')[0]] = kv.split(':')[1]

    df = pd.DataFrame(df)
    for k in mKvExtra:
        df[k]=mKvExtra[k]

    df['libname']=''
    for i in df.index:
        df.ix[i,'libname'] = fxnLibName( df.ix[i] )

    # save to our own key
    df.to_csv(o.libKeyOut, sep='\t', index=False)        

    # save into core template
    if o.coreSheetOut is not None:
        shutil.copyfile( o.coreTemplateIn, o.coreSheetOut )
        ct=openpyxl.load_workbook(filename=o.coreSheetOut)

        wb = ct.active

        # wb['A13']='Sample Name*'

        rowofs=0
        for _,r in df.iterrows():
            shloco=ofsFrom('A18',down=rowofs,right=0)
            wb[ shloco ] = r.libname

            # shloco=ofsFrom('A18',down=rowofs,right=1)
            # wb[ shloco ] = o.coreBuffer

            shloco=ofsFrom('A18',down=rowofs,right=2)
            wb[ shloco ] = float(o.coreConc)

            shloco=ofsFrom('A18',down=rowofs,right=3)
            wb[ shloco ] = float(o.coreVol)

            # shloco=ofsFrom('A18',down=rowofs,right=4)
            # wb[ shloco ] = float(o.coreFraglen)

            shloco=ofsFrom('A18',down=rowofs,right=5)
            bcseq5 = mBcNameSeq_i5i7[ r.p5_barc ][0]
            bcseq7 = mBcNameSeq_i5i7[ r.p7_barc ][1]
            wb[ shloco ]=bcseq7

            shloco=ofsFrom('A18',down=rowofs,right=4)
            wb[ shloco ]=bcseq5

            # shloco=ofsFrom('A18',down=rowofs,right=6)
            # wb[ shloco ] = o.coreShortcode

            shloco=ofsFrom('A18',down=rowofs,right=5)
            wb[ shloco ] = o.coreSpecies

            shloco=ofsFrom('A18',down=rowofs,right=6)
            wb[ shloco ] = "DNA"
            
            if rowofs==0:
                shloco=ofsFrom('A18',down=rowofs,right=7)
                wb[ shloco ] = o.coreNotes


            rowofs+=1

        ct.save(filename=o.coreSheetOut)


if __name__ == '__main__':
    main()
    