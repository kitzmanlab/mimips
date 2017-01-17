cimport cython

from cpython cimport array as c_array
import sys 


# from pysamstats.pyx

DEF BAM_CIGAR_SHIFT=4
DEF BAM_CIGAR_MASK=((1 << BAM_CIGAR_SHIFT) - 1)

DEF BAM_CMATCH     = 0
DEF BAM_CINS       = 1
DEF BAM_CDEL       = 2
DEF BAM_CREF_SKIP  = 3
DEF BAM_CSOFT_CLIP = 4
DEF BAM_CHARD_CLIP = 5
DEF BAM_CPAD       = 6
DEF BAM_CEQUAL     = 7
DEF BAM_CDIFF      = 8


@cython.boundscheck(False) 
@cython.wraparound(False) 
@cython.overflowcheck(False)
cpdef uint32_t clipReadToGapfill_pyx( 
    int gapFillStart, 
    int gapFillEnd, 
    AlignedSegment rd ) except -1:

    cdef bam1_t* src 

    cdef uint32_t cigarNtup 
    cdef uint32_t* pCigar

    cdef cython.str rds  
    cdef c_array.array rdq 

    cdef int cozRefStart, cozRefEnd

    cdef int newCozRefStart

    cdef int cozrefEndClipToGapFill

    cdef uint32_t k, cigop, ciglen
    cdef int32_t nbpungRefSkippedLeft = 0
    cdef int32_t nbpungRefSkippedRight = 0

    cdef bint simple = False

    cdef int32_t nbpungReadPrevBlks = 0 
    cdef int32_t nbpungRefPrevBlks = 0

    # range of indices with respect to reference for gap fill
    cdef int32_t[2] irungrefGapfill
    
    # range of indices with respect to READ for gap fill
    cdef int32_t[2] irungreadGapfill

    cdef list new_cig

    src = rd._delegate
    cigarNtup = pysam_get_n_cigar( src )
    pCigar = pysam_bam_get_cigar(src)

    rds = rd.query_sequence
    rdq = rd.query_qualities

    cozRefStart = rd.reference_start
    cozRefEnd = rd.reference_end-1

    irungreadGapfill[0] = -999999
    irungreadGapfill[1] = -999999

    if (cigarNtup == 0):
        return 0

    elif (cigarNtup == 1) and (pCigar[0] & BAM_CIGAR_MASK) == BAM_CMATCH:
        # easy case: existing cigar str is just a single match block

        #sys.stdout.write('trivial %s:\n'%(rd.qname))
        #sys.stdout.write('  %d cigar %s\n'%(cigarNtup, repr( rd.cigartuples )))

        # rd.ref_end : 0, non-inclusive
        irungrefGapfill[0] = max( cozRefStart, gapFillStart ) - cozRefStart
        irungrefGapfill[1] = min( cozRefEnd, gapFillEnd ) - cozRefStart 

        rd.cigartuples = ( (0,irungrefGapfill[1]-irungrefGapfill[0]+1), )

        #assert irungrefGapfill[0]>=0
        #assert irungrefGapfill[1]>=0

        # debug
        #sys.stdout.write(' 0 irng:%d-%d\n'%(irungrefGapfill[0],irungrefGapfill[1]))

        rd.reference_start = cozRefStart + irungrefGapfill[0]
        rd.query_sequence = rds[ irungrefGapfill[0]:irungrefGapfill[1]+1 ]
        rd.query_qualities = rdq[ irungrefGapfill[0]:irungrefGapfill[1]+1 ]

        return 1

    else:
        # existing cigar string is not a single match block

        new_cig = []

        irungrefGapfill[0] = max( cozRefStart, gapFillStart ) - cozRefStart
        irungrefGapfill[1] = min( cozRefEnd, gapFillEnd ) - cozRefStart

        # debug
        #sys.stdout.write('%s:\n'%(rd.qname))
        #sys.stdout.write('  %d cigar %s\n'%(cigarNtup, repr( rd.cigartuples )))
        #assert irungrefGapfill[0]>=0 and irungrefGapfill[1]>=0

        newCozRefStart = -1

        for k from 0 <= k < cigarNtup:

            # debug
            #sys.stdout.write(' k=%d irungrefGapfill:%d,%d\n'%(k,irungrefGapfill[0],irungrefGapfill[1]))
            #sys.stdout.write(' read prev blks: %d   ref prev blks: %d\n'%(nbpungReadPrevBlks, nbpungRefPrevBlks))

            # if we have seen enough reference bases in previous 
            # blocks to get through gap fill then break
            if nbpungRefPrevBlks - 1 >= irungrefGapfill[1] :
                # debug
                #sys.stdout.write(' breaking\n' )
                break

            cigop = ( pCigar[k] & BAM_CIGAR_MASK )
            ciglen = ( pCigar[k] >> BAM_CIGAR_SHIFT )

            if cigop == BAM_CSOFT_CLIP:
                # don't pass soft-clipped blocks
                nbpungReadPrevBlks += ciglen
            else:

                if (cigop == BAM_CMATCH) or (cigop == BAM_CDEL):
                    # deletion or match - advances on reference

                    #jk debug 12-31-16
                    #if (nbpungRefPrevBlks - 1 < irungrefGapfill[0]) and \
                    #   (nbpungRefPrevBlks + ciglen - 1 >= irungrefGapfill[0])   :
                    if (nbpungRefPrevBlks <= irungrefGapfill[0]) and \
                       (nbpungRefPrevBlks + ciglen - 1 >= irungrefGapfill[0])   :
                        # starts AT OR before gap fill and enters gap fill
 
                        # if current block starts before gap fill then it must 
                        # be match - i.e., we will not start clipped read
                        # with deletion
                        if cigop == BAM_CMATCH:

                            # if this is the first match within the gapfill, record position
                            if newCozRefStart < 0:  # and cigop==BAM_CMATCH
                                newCozRefStart = cozRefStart + irungrefGapfill[0]

                            nbpungRefSkippedLeft = irungrefGapfill[0] - (nbpungRefPrevBlks)

                            irungreadGapfill[0] = nbpungReadPrevBlks + nbpungRefSkippedLeft

                            # debug
                            #sys.stdout.write(' 1start gf l %d skipl=%d\n'%(irungrefGapfill[0],nbpungRefSkippedLeft))

                            # does the current block extend past gap fill?
                            if (nbpungRefPrevBlks + ciglen - 1 > irungrefGapfill[1]):
                                nbpungRefSkippedRight =  nbpungRefPrevBlks + ciglen - irungrefGapfill[1] - 1  # (irungrefGapfill[1] - (nbpungRefPrevBlks + ciglen - 1))

                                irungreadGapfill[1] = nbpungReadPrevBlks + (ciglen - nbpungRefSkippedRight) - 1

                                # debug
                                #sys.stdout.write(' 2past gf r %d skipl=%d skipr=%d\n'%(irungrefGapfill[1],nbpungRefSkippedLeft,nbpungRefSkippedRight))
                            else:
                                nbpungRefSkippedRight = 0                                

                            # debug
                            #sys.stdout.write(' 3appending cigop %d %d\n'%(cigop, ciglen - nbpungRefSkippedLeft - nbpungRefSkippedRight))
                            new_cig.append( (cigop, ciglen - nbpungRefSkippedLeft - nbpungRefSkippedRight) )

                            # if match, increment read position
                            nbpungReadPrevBlks += ciglen 

                    elif (nbpungRefPrevBlks > irungrefGapfill[0]):
                    #jk debug 12-31-16
                    #elif (nbpungRefPrevBlks - 1 >= irungrefGapfill[0]):
                        # starts inside gap fill
                
                        # if this is the first match within the gapfill, record position
                        if newCozRefStart < 0 and  cigop==BAM_CMATCH:
                            newCozRefStart = nbpungRefPrevBlks + cozRefStart

                        # have we not started to pass through blocks yet?
                        #  (eg this block starts exactly at beginning of gap fill, or a deletion block was the first thing to cross into the gf)
                        if irungreadGapfill[0] < 0:
                            irungreadGapfill[0] = nbpungReadPrevBlks

                        # does the current block extend past gap fill?
                        if (nbpungRefPrevBlks + ciglen - 1 > irungrefGapfill[1]):
                            nbpungRefSkippedRight =  nbpungRefPrevBlks + ciglen - irungrefGapfill[1] - 1

                            # debug
                            #sys.stdout.write(' 4past gf r %d, skipr %d \n'%(irungrefGapfill[1],nbpungRefSkippedRight)) 
                        else:
                            nbpungRefSkippedRight = 0

                        # do not end on a deletion
                        if not ( (nbpungRefPrevBlks + ciglen - 1 >= irungrefGapfill[1]) and cigop==BAM_CDEL ):
                            # debug
                            #sys.stdout.write(' 5appending cigop %d %d\n'%(cigop, ciglen - nbpungRefSkippedRight))

                            irungreadGapfill[1] = nbpungReadPrevBlks + (ciglen - nbpungRefSkippedRight) - 1
                            new_cig.append( (cigop, ciglen - nbpungRefSkippedRight) )
                        else:
                            # debug
                            #sys.stdout.write(' skipping since at end; would have been cigop %d %d \n'%(cigop, ciglen - nbpungRefSkippedRight))
                            break

                        if cigop == BAM_CMATCH:
                            # if match, increment read position
                            nbpungReadPrevBlks += ciglen 
                    else:
                        # does not start inside gap fill, does not enter gap fill

                        if cigop == BAM_CMATCH:
                            # if match, increment read position
                            nbpungReadPrevBlks += ciglen 

                    #
                    #
                    # regardless of position wrt gapfill, if it's a match or del, increment ref position
                    nbpungRefPrevBlks += ciglen 
                
                elif cigop == BAM_CINS:

                    if nbpungRefPrevBlks - 1 >= irungrefGapfill[0]:

                        new_cig.append( (cigop, ciglen) )                        

                    # increment read position
                    nbpungReadPrevBlks += ciglen 


        if irungreadGapfill[1] < 0:
            # read did not cross through the end of the gap fill, therefore we never set this
            irungreadGapfill[1] = nbpungReadPrevBlks - 1 


        if newCozRefStart < 0:
            #sys.stdout.write('no match blocks recorded within gf %s\n'%( repr(rd) ))
            return 0

        rd.cigartuples = new_cig

        #MISMATCH_MATE_ALIGNMENT_START irungrefGapfill[0]>=0,irungrefGapfill[0]
        #assert irungreadGapfill[0]>=0,irungreadGapfill[0]
        #assert irungreadGapfill[1]>=0,irungreadGapfill[1]

        rd.query_sequence = rds[ irungreadGapfill[0]:irungreadGapfill[1]+1 ]
        #debugging 12-31-2016
        #rd.reference_start = cozRefStart + irungrefGapfill[0]
        rd.reference_start = newCozRefStart
        rd.query_qualities = rdq[ irungreadGapfill[0]:irungreadGapfill[1]+1 ]

        # debug
        #sys.stdout.write( 'NEW ref_start=%d  cigar=%s\n'%( rd.reference_start, repr( new_cig ) ))

        return 1


