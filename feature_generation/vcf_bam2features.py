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
    parser.add_argument("-r", "--ref-file", action="store", dest="refFile", type=str, required=True, metavar='refFile', help="Full Path to the reference genome file to be used for feature generation")
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
    txt_out = os.path.join(args.outdir,args.outFile)
    #txt_fh = open(txt_out, "wb")
    #txt_fh.write("Tumor_Sample_Barcode\tChromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele1\treads_all\treads_pp\treads_mate_unmapped\treads_mate_other_chr\treads_mate_same_strand\treads_faceaway\treads_softclipped\treads_duplicat\tgc\tmatches\tmismatches\tdeletions\tinsertions\tA/C/T/G/N\tmean_tlen\trms_tlen\tstd_tlen\tread_mapq0\trms_mapq\tmax_mapq\trms_baseq\trms_baseq_matches\trms_baseq_mismatches\n")
    rec_dict_list = []
    keys = []
    keyorder = ['chrom','pos','ref','reads_all','reads_fwd','reads_rev','reads_pp','reads_pp_fwd','reads_pp_rev','matches','matches_fwd','matches_rev','matches_pp','matches_pp_fwd','matches_pp_rev','mismatches','mismatches_fwd','mismatches_rev','mismatches_pp','mismatches_pp_fwd','mismatches_pp_rev','deletions','deletions_fwd','deletions_rev','deletions_pp','deletions_pp_fwd','deletions_pp_rev','insertions','insertions_fwd','insertions_rev','insertions_pp','insertions_pp_fwd','insertions_pp_rev','A','A_fwd','A_rev','A_pp','A_pp_fwd','A_pp_rev','C','C_fwd','C_rev','C_pp','C_pp_fwd','C_pp_rev','G','G_fwd','G_rev','G_pp','G_pp_fwd','G_pp_rev','T','T_fwd','T_rev','T_pp','T_pp_fwd','T_pp_rev','N','N_fwd','N_rev','N_pp','N_pp_fwd','N_pp_rev']
# iterate over statistics, one record at a time
    for record in vcf_reader:
        chromosome = record.CHROM
        position = int(record.POS)
        ref = record.REF
        alt = record.ALT[0]
        for rec in pysamstats.stat_variation_strand(bam_to_process, args.refFile, chrom=chromosome, start=position, end=position+1,truncate=True):
            rec = collections.OrderedDict(sorted(rec.items(),key=lambda i:keyorder.index(i[0])))
            rec = MyOrderedDict(rec)
            #rec.prepend('Tumor_Seq_Allele1',alt)
            #rec.prepend('Reference_Allele',ref)
            #rec.prepend('Start_Position',position)
            rec.prepend('Tumor_Sample_Barcode',args.sampleName)
            #print rec
            keys=rec.keys()
            print "Org:",chromosome,position,ref,alt,rec['chrom'],rec['pos'],rec['ref'],"\n"
            rec_dict_list.append(rec)
    #Write output
    with open(txt_out, 'wb') as output_file:
        dict_writer = csv.DictWriter(output_file, keys, delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(rec_dict_list)
    return
    
class MyOrderedDict(collections.OrderedDict):
    def prepend(self, key, value, dict_setitem=dict.__setitem__):
        root = self._OrderedDict__root
        first = root[1]
        if key in self:
            link = self._OrderedDict__map[key]
            link_prev, link_next, _ = link
            link_prev[1] = link_next
            link_next[0] = link_prev
            link[0] = root
            link[1] = first
            root[1] = first[0] = link
        else:
            root[1] = first[0] = self._OrderedDict__map[key] = [root, first, key]
            dict_setitem(self, key, value)
            
# Run the whole script
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("vcf_bam2features: Elapsed time was %g seconds", totaltime)
    sys.exit(0)