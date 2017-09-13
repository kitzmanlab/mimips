from libc.stdint cimport int32_t, uint32_t, uint8_t, uint64_t, int64_t
from pysam.libchtslib cimport bam1_t #, bam_pileup1_t
from pysam.libcalignmentfile cimport AlignmentFile #, IteratorRowRegion
from pysam.libcalignedsegment cimport  AlignedSegment, pysam_get_n_cigar, pysam_bam_get_cigar
