from setuptools import setup
#from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import pysam
import numpy 
import glob
import os.path as op

setup(
	name="mimips",

    # cmdclass = {'build_ext':build_ext},

    packages=['mimips','mimips.uniformity'],

    ext_modules = cythonize(Extension("mip_pipe_core", ["mimips/mip_pipe_core.pyx"]),
              annotate=True),

    include_dirs = [numpy.get_include()]+pysam.get_include(),

    entry_points = {
    	'console_scripts': [ 'rev_read_strip_mid = mimips.rev_read_strip_mid:main',
    		 	 	         'merged_pe_strip_mid = mimips.merged_pe_strip_mid:main',
                            'est_libsize_already_dedupd = mimips.est_libsize_already_dedupd:main',
    		 	 	         'flag_dup_mid_reads = mimips.flag_dup_mid_reads:main',
    		 	 	         'exact_trim_mip_arms = mimips.exact_trim_mip_arms_v2:main',
    		 	 	         'annotate_bam_by_mip_SE = mimips.annotate_bam_by_mip_SE:main',
    		 	 	         'annotate_bam_by_mip = mimips.annotate_bam_by_mip:main',
    		 	 	         'gather_target_coverage_at_thresholds = mimips.uniformity.gather_target_coverage_at_thresholds:main',
                             'gather_perbp_uniformity = mimips.uniformity.gather_perbp_uniformity:main',
    		 	 	         'probe_hit_count_from_bam = mimips.uniformity.probe_hit_count_from_bam:main',
    		 	 	         'target_coverage_uniformity = mimips.uniformity.target_coverage_uniformity:main',
                             'forcore_mip_sample_sheet = mimips.forcore_mip_sample_sheet:main',
                             'join_mip_libsheet_core_iemsheet = mimips.join_mip_libsheet_core_iemsheet:main',
                             'sex_chrom_check = mimips.uniformity.sex_chrom_checks:main'
    		 	 	          ]
    }
)





