### vcf_bam2features
```
usage: vcf_bam2features.py [options]

This module gets the detailed bam information for a given chromosome, start,
end, reference allele and alternate allele and 1bp flanking from a vcf

optional arguments:
  -h, --help            show this help message and exit
  -o OutFile, --output-file OutFile
                        Name of the output file
  -b BamFile, --bam-file BamFile
                        Full Path to the bam file to be used for feature
                        generation
  -v, --verbose         make lots of noise
  -d /somepath/output, --output-dir /somepath/output
                        Full Path to the output directory.
  -s SomeID, --sample-name SomeID
                        Sample ID for the bam file, if not provided the prefix
                        of the bam file will be used
  -i /somepath/variants.vcf, --input-vcf /somepath/variants.vcf
                        Full Path to variant vcf
```