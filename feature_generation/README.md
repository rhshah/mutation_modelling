### vcf_bam2features

- This module gets the detailed bam information for a given chromosome, start, end, reference allele and alternate allele
- In test for a vcf with ~300k variants using 40 processors it takes ~5 min for 1 bam file.

#### Requirements:
- pyvcf : [v0.6.7](http://pyvcf.readthedocs.io/en/latest/INTRO.html)
- pandas : [v0.20.3](http://pandas.pydata.org/)
- pysam : [v0.11.2.2](https://github.com/pysam-developers/pysam)
- pysamstats : [v1.0.1](https://github.com/alimanfoo/pysamstats)
- joblib : [v0.11](https://pythonhosted.org/joblib/)
- coloredlogs: [v7.1](https://coloredlogs.readthedocs.io/en/latest/)

## Install using virtualenv
Using [virtualenv](https://virtualenv.pypa.io) you can do the following:

```
virtualenv -p /path/to/python/python-2.7.10/bin/python mutation_modelling
cd mutation_modelling
source bin/activate
pip install coloredlogs==7.1
pip install joblib==0.11
pip install pandas==0.20.3
pip install pyvcf==0.6.7
pip install pysam==0.11.2.2
pip install pysamstats==1.0.1
git clone https://github.com/rhshah/mutation_modelling.git
python mutation_modelling/feature_generation/vcf_bam2features.py -h

```

#### Usage
```
usage: vcf_bam2features.py [options]

This module gets the detailed bam information for a given chromosome, start,
end, reference allele and alternate allele

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

#### Example

```
python vcf_bam2features.py -b /path/to/some.bam -v -d /home/shahr2/virtualenv/mutation_modelling/ -s AX279423-T -i /path/to/some.vcf -o test1 -r /path/to/genome.fasta -p 20 -mapq 0 -baseq 0 
```

#### Output
It will produce 5 files.
1. prefix_variant-stats.txt
2. prefix_baseq-stats.txt
3. prefix_mapq-stats.txt
4. prefix_gc-stats.txt
5. prefix_merged-stats.txt which is essentially combination of all

Note: For details of each column please refer the [pysamstats](https://github.com/alimanfoo/pysamstats) documentation.