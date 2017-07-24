"""
vcf_bam2features
~~~~~~~~~~~~~~~~~~~~~

:Description: This module gets the detailed bam information for a given chromosome, start, end, reference allele and alternate allele and 1bp flanking from a vcf

"""

'''
Created on July 26, 2017
Description: This module gets the detailed bam information for a given chromosome, start, end, reference allele and alternate allele and 1bp flanking from a vcf
@author: Ronak H Shah
::Input::
vcf: location of vcf file for variant sites
bam: bam file from which features are to be generated
ouput_dir: out directory if not given current working directory
output_file: output file
::Output::


::Example Run::
```
python vcf_bam2features.py 
```

'''
import argparse
import time
import os
import sys
import logging
import glob
import textwrap

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('vcf_bam2features')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("vcf_bam2features: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass

try:
    import pysam
except ImportError:
    logger.fatal("vcf_bam2features: pysam is not installed, please install pysam as it is required to generate features. pysam version == 0.8.4")
    sys.exit(1)
try:
    import pysamstats
except ImportError:
    logger.fatal("vcf_bam2features: pysamstats is not installed, please install pysamstats as it is required to generate feature. pysamstats version == 0.24.3")
    sys.exit(1)


try:
    import vcf
except ImportError:
    logger.fatal("vcf_bam2features: pyvcf is not installed, please install pyvcf as it is required to run the mapping. pyvcf version == 0.6.7")
    sys.exit(1)


#Run all sub function    
def main():
    parser = argparse.ArgumentParser(prog='vcf_bam2features.py', description='This module gets the detailed bam information for a given chromosome, start, end, reference allele and alternate allele and 1bp flanking from a vcf', usage='%(prog)s [options]')
    parser.add_argument("-o", "--output-file", action="store", dest="outFile", type=str, required=True, metavar='OutFile', help="Name of the output file")
    parser.add_argument("-b", "--bam-file", action="store", dest="bamFile", type=str, required=True, metavar='BamFile', help="Full Path to the bam file to be used for feature generation")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
    parser.add_argument("-d", "--output-dir", action="store", dest="outdir", type=str, required=False, metavar='/somepath/output', help="Full Path to the output directory.")
    parser.add_argument("-s", "--sample-name", action="store", dest="sampleName", type=str, required=False, metavar='SomeID', help="Sample ID for the bam file, if not provided the prefix of the bam file will be used")
    parser.add_argument("-i", "--input-vcf", action="store", dest="inputVcf", type=str, required=False,metavar='/somepath/variants.vcf', help="Full Path to variant vcf")
    args = parser.parse_args()
	   
    args = validate_inputs(args)
    
    generate_features(args)
    if(args.verbose):
        logger.info("vcf_bam2features: Finished generating flanking sequence")
        #logger.info("vcf_bam2features: Flanking sequence is written in %s", outfile)


#Validate Given Input Arguments
def validate_inputs(args):
    if(os.path.isfile(args.inputVcf)):
        if(args.verbose):
            logger.info("vcf_bam2features: Using %s as the vcf file for selection of sites", args.inputVcf)
        pass
    else:
        if(args.verbose):
            logger.critical("vcf_bam2features: %s is not a valid file please rerun using a valid vcf file", args.inputVcf)
        sys.exit(1)
    
    if(os.path.isfile(args.bamFile)):
        if(args.verbose):
            logger.info("vcf_bam2features: Using %s as the bam file to generate features", args.bamFile)
        pass
    else:
        if(args.verbose):
            logger.critical("vcf_bam2features: %s is not a valid file please rerun using a valid bam file", args.bamFile)
        sys.exit(1)
    if(args.outdir):
        if(os.path.isdir(args.outdir)):
            if(args.verbose):
                logger.info("vcf_bam2features: %s is a valid directory and will be used to write output", args.outdir)
            pass
        else:
            if(args.verbose):
                logger.critical("vcf_bam2features: %s is not a valid directory please rerun using a valid directory", args.outdir)
            sys.exit(1)
    else:
        args.outdir = os.getcwd()
        logger.info("vcf_bam2features: Output directory is not specified and thus %s will be used to write output", args.outdir)
    
    if(type(args.outFile) is str):
        if(args.verbose):
            logger.info("vcf_bam2features: %s is a string and will be used as the output file name", args.outFile)
        pass
    else:
        logger.info("vcf_bam2features: %s is not a string please use a proper string as a name", args.outFile)
        sys.exit(1)
    
    if(args.sampleName):
        if(type(args.sampleName) is str):
            if(args.verbose):
                logger.info("vcf_bam2features: %s is a string and will be used in output file", args.sampleName)
        pass
    else:
        args.sampleName = os.path.splitext(args.bamFile)[0]
        if(args.verbose):
            logger.info("vcf_bam2features: %s is a string and will be used in output file", args.sampleName)
    
    return(args)
def generate_features(args):
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))
    bam_to_process = pysam.AlignmentFile(args.bamFile)
# iterate over statistics, one record at a time
    for record in vcf_reader:
        chromosome = record.CHROM
        position = int(record.POS)
        if(len(str(record.REF)) > len(str(record.ALT))):
            start = position
            end = position + 2
        else:
            start = position -1
            end = position + 1
        for rec in pysamstats.stat_coverage(bam_to_process, chrom=chromosome, start=start, end=end):
            print 
            args.sampleName, chromosome, record.REF, record.ALT,
            rec['chrom'], rec['pos'], 
            rec['reads_all'], rec['reads_all_fwd'], rec['reads_all_rev'], 
            rec['reads_pp'], rec['reads_pp_fwd'], rec['reads_pp_rev'], 
            rec['reads_mate_unmapped'], rec['reads_mate_unmapped_fwd'], rec['reads_mate_unmapped_rev'], 
            rec['reads_mate_other_chr'], rec['reads_mate_other_chr_fwd'], rec['reads_mate_other_chr_rev'],
            rec['reads_mate_same_strand'], rec['reads_mate_same_strand_fwd'], rec['reads_mate_same_strand_rev'],
            rec['reads_faceaway'], rec['reads_faceaway_fwd'], rec['reads_faceaway_rev'],
            rec['reads_softclipped'], rec['reads_softclipped_fwd'], rec['reads_softclipped_rev'],
            rec['reads_duplicate'], rec['reads_duplicate_fwd'], rec['reads_duplicate_rev'],
            rec['gc'], 
            rec['matches'], rec['matches_pp'], 
            rec['mismatches'], rec['mismatches_pp'], 
            rec['deletions'], rec['deletions_pp'], 
            rec['insertions'], rec['insertions_pp'], 
            rec['A'], rec['A_pp'], 
            rec['C'], rec['C_pp'], 
            rec['T'], rec['T_pp'], 
            rec['G'], rec['G_pp'], 
            rec['N'], rec['N_pp'],
            rec['mean_tlen'], rec['mean_tlen_pp'], rec['rms_tlen'], rec['rms_tlen_pp'], rec['std_tlen'], rec['std_tlen_pp'],  
            rec['reads_mapq0'], rec['rms_mapq'], rec['rms_mapq_pp'], rec['max_mapq'], rec['max_mapq_pp'], 
            rec['rms_baseq'], rec['rms_baseq_matches'], rec['rms_baseq_matches_pp'], rec['rms_baseq_mismatches'], rec['rms_baseq_mismatches_pp']
            
    return()
# Run the whole script
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("vcf_bam2features: Elapsed time was %g seconds", totaltime)
    sys.exit(0)