import sys
import argparse
from collections import defaultdict

import pysam

import numpy as np

def main():

    opts = argparse.ArgumentParser()

    opts.add_argument('--inBam',default=None,dest='inBam')
    opts.add_argument('--outBam',default=None,dest='outBam')

    opts.add_argument('--sampleN',type=int,default=100,dest='sampleN')
    
    o = opts.parse_args()

    ######################################################

    # strategy - step thourgh and add lines to queue
    # everytime a gap of >= MAX_GAP between consecutive start 
    # positions in encountered, process the queued lines,
    # sampling by probe id 

    MAX_GAP = 10

    bamIn = pysam.Samfile( o.inBam, 'rb' )
    bamOut = pysam.Samfile( o.outBam, 'wb', template=bamIn )

    lreads = []
    m_pbidx_li_reads = defaultdict(list)

    last_pos = -10000
    last_refid = None

    for r in bamIn:

        if r.reference_start >= last_pos + MAX_GAP - 1 or r.reference_id != last_refid:

            # sys.stderr.write('resetting queue @ {}\n'.format(r.reference_start))

            # clear old queue
            lout = []
            for pbidx in m_pbidx_li_reads:
                # are there more reads for this probe than our max per-probe sampling depth?
                if len(m_pbidx_li_reads[pbidx]) > o.sampleN:
                    # sample them
                    li_chosen = np.random.choice( 
                        m_pbidx_li_reads[pbidx], 
                        o.sampleN,
                        replace=False )
                    lout += [ lreads[j] for j in li_chosen ]
                else:
                    # add all reads
                    lout += [ lreads[j] for j in m_pbidx_li_reads[pbidx] ]

            for l in sorted(lout,key=lambda ol:ol.reference_start): 
                bamOut.write(l)

            # clear queue
            lreads = []
            m_pbidx_li_reads = defaultdict(list)


        last_pos = r.reference_start
        last_refid = r.reference_id

        lreads.append( r )
        m_pbidx_li_reads[ dict(r.tags)['ZP'] ].append( len(lreads)-1 )


    if len(lreads)>0:
        # clear old queue
        lout = []
        for pbidx in m_pbidx_li_reads:
            # are there more reads for this probe than our max per-probe sampling depth?
            if len(m_pbidx_li_reads[pbidx]) > o.sampleN:
                # sample them
                li_chosen = np.random.choice( 
                    m_pbidx_li_reads[pbidx], 
                    o.sampleN,
                    replace=False )
                lout += [ lreads[j] for j in li_chosen ]
            else:
                # add all reads
                lout += [ lreads[j] for j in m_pbidx_li_reads[pbidx] ]

        for l in sorted(lout,key=lambda ol:ol.reference_start): 
            bamOut.write(l)



if __name__ == '__main__':                
    main()