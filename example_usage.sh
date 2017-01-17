# MIMIPS_PIPELINE_FILES should be an environment variable pointing to the location of the
# mimips "pipeline_files" directory


# 1. For submission to U-M core: generate a sample key (excel formatted) as well as a TSV sample key for downstream processing.

# relevant arguments:
    # --excelSheetIn [path to excel sheet exported from google drive ]
    # --libKeyOut [path to a new TSV file which will go into the current pipeline]
    # --coreTemplateIn [path to blank U-M core sample subsmission sheet]
    # 

forcore_mip_sample_sheet \
    --excelSheetIn capture_plate2.xlsx \
    --libKeyOut key_for_pipeline.tsv \
    --barcodeFasta  ${MIMIPS_PIPELINE_FILES}/kitzman_barcodes.fa \
    --coreTemplateIn IlluminaSampleClient.xlsx \
    --coreSheetOut ~/Google\ Drive/KitzmanLabShared/Shared\ Notebooks/Aortic\ Aneurysm\ MIPS/Core\ submission\ notes/kitzman-CVC-35.xlsx \
    --extraCols "captureset:plate1"\
    --coreConc XXXXXXX \
    --coreVol 40 --coreFraglen 300 --coreShortcode XXXXXX --coreSpecies human --coreBuffer TE\
    --fxnLibName 'lambda r:"%s_%s_%s_%s"%( r.captureset, r.well, r.source_plate, r.source_well ) ' \
    --coreNote "This is a pooled MIPS amplicon library.  Please sequence with paired index reads (10 bp each) and provide index reads as separate FASTQs. Please use enclosed forward sequencing primer (kitzman-MIP-FWD CAGGACGTCAGATGTTATCGAGGTCCGAC), enclosed index1 primer (kitzman-MIP-IDX GTTGGAGGCTCATCGTTCCTCGATACGG) and enclosed read 2 primer (kitzman-MIP-REV CCGTATCGAGGAACGATGAGCCTCCAAC)"


# 2. after data return, pair up received fastq files with the expected libraries and generate a library key with
#    fastq paths for the pipeline

# --libKeyIn [path to TSV file from step 1]
# --libKeyOut [path to new TSV file to make]
# --coreIEMSheet [path to CSV sheet included with illumina data from core]
# --baseDir [path under which directories w/ fastqs can be found]

join_mip_libsheet_core_iemsheet \
    --libKeyIn key_for_pipeline.tsv \
    --libKeyOut key_for_pipeline_withpaths.tsv \
    --coreIEMSheet /path/to/Run_NNNN_piname.csv \
    --baseDir /path/to/piname/


# 3. move to a new directory for processing

pushd /path/to/processing_dir

# 3a. make a copy of snakemake processing pipeline
cp ${MIMIPS_PIPELINE_FILES}/config.yaml ./
cp ${MIMIPS_PIPELINE_FILES}/preprocess_and_align.snake ./


# 3b. modify the yaml file to point to your specific paths (bwa index location, 
#     and your own table of probe locations )

# 3c. run snakemake, pointing it to the sample key we just made, 
#     specifying to run in the current directory, 

source activate mipspy3
snakemake -s preprocess_and_align.snake \
    --configfile config.yaml \
    --config sample_table=key_for_pipeline_withpaths.tsv \
    work_dir=`pwd` \
    prefix=mymipset \
     --cores 16 \
    --printshellcmds 






# TODO FIX acap BUG
# TODO POST to a github