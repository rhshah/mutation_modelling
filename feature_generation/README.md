### vcf_bam2features

#### Requirements:
- pyvcf : [v0.6.7](http://pyvcf.readthedocs.io/en/latest/INTRO.html)
- pandas : [v0.20.3](http://pandas.pydata.org/)
- pysam : [v0.11.2.2](https://github.com/pysam-developers/pysam)
- pysamstats : [v1.0.1](https://github.com/alimanfoo/pysamstats)
- joblib : [v0.11](https://pythonhosted.org/joblib/)
- coloredlogs: [v7.1](https://coloredlogs.readthedocs.io/en/latest/)

#### Usage
```
usage: vcf_bam2features.py [options]

This module gets the detailed bam information for a given chromosome, start,
end, reference allele and alternate allele and 1bp flanking from a vcf

optional arguments:
  -h, --help            show this help message and exit
  -o OutFile, --output-file-prefix OutFile
                        Prefix of the output file
  -b BamFile, --bam-file BamFile
                        Full Path to the bam file to be used for feature
                        generation
  -r refFile, --ref-file refFile
                        Full Path to the reference genome file to be used for
                        feature generation
  -v, --verbose         make lots of noise
  -d /somepath/output, --output-dir /somepath/output
                        Full Path to the output directory.
  -s SomeID, --sample-name SomeID
                        Sample ID for the bam file, if not provided the prefix
                        of the bam file will be used
  -i /somepath/variants.vcf, --input-vcf /somepath/variants.vcf
                        Full Path to variant vcf
  -p 10, --processors 10
                        Total number of processors to use
  -mapq 20, --min-mapq 20
                        Minimum Mapping Quality
  -baseq 20, --min-baseq 20
                        Minimum Base Quality
```