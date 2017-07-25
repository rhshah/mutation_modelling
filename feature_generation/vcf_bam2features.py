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
import csv
import collections
import pprint

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
    import pandas as pd
except ImportError:
    logger.fatal("vcf_bam2features: pandas is not installed, please install pandas as it is required to generate features. pandas version == 0.20.3")
    sys.exit(1)
try:
    from joblib import Parallel, delayed
except ImportError:
    logger.fatal("vcf_bam2features: joblib is not installed, please install joblib as it is required to generate features. joblib version == 0.11")
    sys.exit(1)
try:
    import vcf
except ImportError:
    logger.fatal("vcf_bam2features: pyvcf is not installed, please install pyvcf as it is required to run the mapping. pyvcf version == 0.6.7")
    sys.exit(1)


#Run all sub function    
def main():
    parser = argparse.ArgumentParser(prog='vcf_bam2features.py', description='This module gets the detailed bam information for a given chromosome, start, end, reference allele and alternate allele and 1bp flanking from a vcf', usage='%(prog)s [options]')
    parser.add_argument("-o", "--output-file-prefix", action="store", dest="outFile", type=str, required=True, metavar='OutFile', help="Prefix of the output file")
    parser.add_argument("-b", "--bam-file", action="store", dest="bamFile", type=str, required=True, metavar='BamFile', help="Full Path to the bam file to be used for feature generation")
    parser.add_argument("-r", "--ref-file", action="store", dest="refFile", type=str, required=True, metavar='refFile', help="Full Path to the reference genome file to be used for feature generation")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
    parser.add_argument("-d", "--output-dir", action="store", dest="outdir", type=str, required=False, metavar='/somepath/output', help="Full Path to the output directory.")
    parser.add_argument("-s", "--sample-name", action="store", dest="sampleName", type=str, required=True, metavar='SomeID', help="Sample ID for the bam file, if not provided the prefix of the bam file will be used")
    parser.add_argument("-i", "--input-vcf", action="store", dest="inputVcf", type=str, required=True,metavar='/somepath/variants.vcf', help="Full Path to variant vcf")
    parser.add_argument("-p", "--processors", action="store", dest="processors", type=int, default=10, required=False, metavar='10', help="Total number of processors to use")
    
    args = parser.parse_args()
	   
    args = validate_inputs(args)
    
    generate_features(args.inputVcf,args.sampleName,args.bamFile,args.refFile,args.outdir,args.outFile,args.processors)
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
    if(os.path.isfile(args.refFile)):
        if(args.verbose):
            logger.info("vcf_bam2features: Using %s as the reference file to generate features", args.bamFile)
        pass
    else:
        if(args.verbose):
            logger.critical("vcf_bam2features: %s is not a valid file please rerun using a valid reference file", args.refFile)
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
            logger.info("vcf_bam2features: %s is a string and will be used as the output file name prefix", args.outFile)
        pass
    else:
        logger.info("vcf_bam2features: %s is not a string please use a proper string as a name prefix", args.outFile)
        sys.exit(1)
    if(type(args.processors) is int):
        if(args.verbose):
            logger.info("vcf_bam2features: %s is a int and will be used as number of processors", args.processors)
        pass
    else:
        logger.info("vcf_bam2features: %s is not a int please use a proper integer values for processors", args.processors)
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

#Generate Features
def generate_features(inputVcf,sampleName,bamFile,refFile,outdir,outFile,processors):
    pp = pprint.PrettyPrinter(indent=4)
    keyorder_a = ['Tumor_Sample_Barcode','chrom','pos','ref','alt','reads_all','reads_fwd','reads_rev','reads_pp','reads_pp_fwd','reads_pp_rev','matches','matches_fwd','matches_rev','matches_pp','matches_pp_fwd','matches_pp_rev','mismatches','mismatches_fwd','mismatches_rev','mismatches_pp','mismatches_pp_fwd','mismatches_pp_rev','deletions','deletions_fwd','deletions_rev','deletions_pp','deletions_pp_fwd','deletions_pp_rev','insertions','insertions_fwd','insertions_rev','insertions_pp','insertions_pp_fwd','insertions_pp_rev','A','A_fwd','A_rev','A_pp','A_pp_fwd','A_pp_rev','C','C_fwd','C_rev','C_pp','C_pp_fwd','C_pp_rev','G','G_fwd','G_rev','G_pp','G_pp_fwd','G_pp_rev','T','T_fwd','T_rev','T_pp','T_pp_fwd','T_pp_rev','N','N_fwd','N_rev','N_pp','N_pp_fwd','N_pp_rev']
    keyorder_b = ['Tumor_Sample_Barcode','chrom','pos','ref','alt','reads_all','reads_fwd','reads_rev','reads_pp','reads_pp_fwd','reads_pp_rev','matches','matches_fwd','matches_rev','matches_pp','matches_pp_fwd','matches_pp_rev','mismatches','mismatches_fwd','mismatches_rev','mismatches_pp','mismatches_pp_fwd','mismatches_pp_rev','rms_baseq','rms_baseq_fwd','rms_baseq_rev','rms_baseq_pp','rms_baseq_pp_fwd','rms_baseq_pp_rev','rms_baseq_matches','rms_baseq_matches_fwd','rms_baseq_matches_rev','rms_baseq_matches_pp','rms_baseq_matches_pp_fwd','rms_baseq_matches_pp_rev','rms_baseq_mismatches','rms_baseq_mismatches_fwd','rms_baseq_mismatches_rev','rms_baseq_mismatches_pp','rms_baseq_mismatches_pp_fwd','rms_baseq_mismatches_pp_rev']
    keyorder_c = ['Tumor_Sample_Barcode','chrom','pos','reads_all','reads_fwd','reads_rev','reads_pp','reads_pp_fwd','reads_pp_rev','reads_mapq0','reads_mapq0_fwd','reads_mapq0_rev','reads_mapq0_pp','reads_mapq0_pp_fwd','reads_mapq0_pp_rev','rms_mapq','rms_mapq_fwd','rms_mapq_rev','rms_mapq_pp','rms_mapq_pp_fwd','rms_mapq_pp_rev','max_mapq','max_mapq_fwd','max_mapq_rev','max_mapq_pp','max_mapq_pp_fwd','max_mapq_pp_rev']
    vcf_reader = vcf.Reader(open(inputVcf, 'r'))
    txt_out1 = os.path.join(outdir,outFile,"_variation.txt")
    txt_out2 = os.path.join(outdir,outFile,"_baseq.txt")
    txt_out3 = os.path.join(outdir,outFile,"_mapq.txt")
    #rec_variation_df_list = []
    #iterate over statistics, one record at a time
    rec_variation_df_list = Parallel(n_jobs=processors, backend="threading")(delayed(run_pysamstats_variation)(bamFile,refFile,sampleName,record)
                           for record in vcf_reader)
    #print "typeof",type(rec_variation_df_list),"\n"
    pp.pprint(rec_variation_df_list[0:10])
    logger.info("Total Record in list of df:%s", len(rec_variation_df_list))
    df1 = pd.DataFrame.from_dict(rec_variation_df_list)
    df1.to_csv(txt_out1,sep="\t",ignore_index=True)
    
    #rec_baseq_df_list = []
    #iterate over statistics, one record at a time
    rec_baseq_df_list = Parallel(n_jobs=processors, backend="threading")(delayed(run_pysamstats_baseq)(bamFile,refFile,sampleName,record)
                           for record in vcf_reader)
    #print "typeof",type(rec_baseq_df_list),"\n"
    logger.info("Total Record in list of df:%s", len(rec_baseq_df_list))
    df2 = pd.DataFarame.from_dict(rec_baseq_df_list)
    df2.to_csv(txt_out2,sep="\t",ignore_index=True)
    
    #rec_mapq_df_list = []
    #iterate over statistics, one record at a time
    rec_mapq_df_list = Parallel(n_jobs=processors, backend="threading")(delayed(run_pysamstats_mapq)(bamFile,refFile,sampleName,record)
                           for record in vcf_reader)
    #print "typeof",type(rec_mapq_df_list),"\n"
    logger.info("Total Record in list of df:%s", len(rec_mapq_df_list))
    df3 = pd.DataFrame(rec_mapq_df_list)
    df3.to_csv(txt_out3,sep="\t",ignore_index=True)
    #with open(txt_out, 'wb') as output_file:
    #    dict_writer = csv.DictWriter(output_file, keys, delimiter='\t')
    #    dict_writer.writeheader()
    #    dict_writer.writerows(rec_dict_list)
    return

#Run PySamStats
def run_pysamstats_variation(bamFile,refFile,sampleName,record):
    bam_to_process = pysam.AlignmentFile(bamFile)
    keyorder = ['Tumor_Sample_Barcode','chrom','pos','ref','alt','reads_all','reads_fwd','reads_rev','reads_pp','reads_pp_fwd','reads_pp_rev','matches','matches_fwd','matches_rev','matches_pp','matches_pp_fwd','matches_pp_rev','mismatches','mismatches_fwd','mismatches_rev','mismatches_pp','mismatches_pp_fwd','mismatches_pp_rev','deletions','deletions_fwd','deletions_rev','deletions_pp','deletions_pp_fwd','deletions_pp_rev','insertions','insertions_fwd','insertions_rev','insertions_pp','insertions_pp_fwd','insertions_pp_rev','A','A_fwd','A_rev','A_pp','A_pp_fwd','A_pp_rev','C','C_fwd','C_rev','C_pp','C_pp_fwd','C_pp_rev','G','G_fwd','G_rev','G_pp','G_pp_fwd','G_pp_rev','T','T_fwd','T_rev','T_pp','T_pp_fwd','T_pp_rev','N','N_fwd','N_rev','N_pp','N_pp_fwd','N_pp_rev']  
    chromosome = record.CHROM
    position = record.POS
    ref = record.REF
    alt = record.ALT[0]
    if(len(str(ref)) > len(str(alt))):
        start = position
        end = position + 1
    else:
        start = position - 1
        end = position 
    for rec in pysamstats.stat_variation_strand(bam_to_process, refFile, chrom=chromosome, start=start, end=end, one_based=True, truncate=True):
        rec['alt'] = alt
        rec['pos']=position
        rec['Tumor_Sample_Barcode']=sampleName
        rec = collections.OrderedDict(sorted(rec.items(),key=lambda i:keyorder.index(i[0])))
        #print "Org:",chromosome,position,ref,alt,rec['chrom'],rec['pos'],rec['ref'],"\n"
        return(rec)
def run_pysamstats_baseq(bamFile,refFile,sampleName,record):
    bam_to_process = pysam.AlignmentFile(bamFile)
    keyorder = ['Tumor_Sample_Barcode','chrom','pos','ref','alt','reads_all','reads_fwd','reads_rev','reads_pp','reads_pp_fwd','reads_pp_rev','matches','matches_fwd','matches_rev','matches_pp','matches_pp_fwd','matches_pp_rev','mismatches','mismatches_fwd','mismatches_rev','mismatches_pp','mismatches_pp_fwd','mismatches_pp_rev','rms_baseq','rms_baseq_fwd','rms_baseq_rev','rms_baseq_pp','rms_baseq_pp_fwd','rms_baseq_pp_rev','rms_baseq_matches','rms_baseq_matches_fwd','rms_baseq_matches_rev','rms_baseq_matches_pp','rms_baseq_matches_pp_fwd','rms_baseq_matches_pp_rev','rms_baseq_mismatches','rms_baseq_mismatches_fwd','rms_baseq_mismatches_rev','rms_baseq_mismatches_pp','rms_baseq_mismatches_pp_fwd','rms_baseq_mismatches_pp_rev']
    chromosome = record.CHROM
    position = record.POS
    ref = record.REF
    alt = record.ALT[0]
    if(len(str(ref)) > len(str(alt))):
        start = position
        end = position + 1
    else:
        start = position - 1
        end = position 
    for rec in pysamstats.stat_baseq_ext(bam_to_process, refFile, chrom=chromosome, start=start, end=end,one_based=True,truncate=True):
        rec['alt'] = alt
        rec['pos']=position
        rec['Tumor_Sample_Barcode']=sampleName
        rec = collections.OrderedDict(sorted(rec.items(),key=lambda i:keyorder.index(i[0])))
        #print "Org:",chromosome,position,ref,alt,rec['chrom'],rec['pos'],rec['ref'],"\n"
        return(rec)
def run_pysamstats_mapq(bamFile,refFile,sampleName,record):
    bam_to_process = pysam.AlignmentFile(bamFile)
    keyorder = ['Tumor_Sample_Barcode','chrom','pos','reads_all','reads_fwd','reads_rev','reads_pp','reads_pp_fwd','reads_pp_rev','reads_mapq0','reads_mapq0_fwd','reads_mapq0_rev','reads_mapq0_pp','reads_mapq0_pp_fwd','reads_mapq0_pp_rev','rms_mapq','rms_mapq_fwd','rms_mapq_rev','rms_mapq_pp','rms_mapq_pp_fwd','rms_mapq_pp_rev','max_mapq','max_mapq_fwd','max_mapq_rev','max_mapq_pp','max_mapq_pp_fwd','max_mapq_pp_rev']
    chromosome = record.CHROM
    position = record.POS
    ref = record.REF
    alt = record.ALT[0]
    if(len(str(ref)) > len(str(alt))):
        start = position
        end = position + 1
    else:
        start = position - 1
        end = position 
    for rec in pysamstats.stat_mapq_strand(bam_to_process, refFile, chrom=chromosome, start=start, end=end,one_based=True,truncate=True):
        rec['alt'] = alt
        rec['pos']=position
        rec['Tumor_Sample_Barcode']=sampleName
        rec = collections.OrderedDict(sorted(rec.items(),key=lambda i:keyorder.index(i[0])))
        #print "Org:",chromosome,position,ref,alt,rec['chrom'],rec['pos'],rec['ref'],"\n"
        return(rec)

# Run the whole script
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("vcf_bam2features: Elapsed time was %g seconds", totaltime)
    sys.exit(0)