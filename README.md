# mergekit
NGS Bioinformatics / Toolkit to assist overlapping paired-end reads merging


## Scripts

### check_merge.pl

**Development script not for production**

### merge.sh

Script that produces the merged files starting from the FASTQ pairs in the `datasets` directory. 
Requires merging tools (namely _flash_, _NGmerge_, _vsearch_ and optionally _usearch_)
to be available, or will attempt to install them via conda.

### detect_region.pl

Will check an input sequence in FASTA or FASTQ format and detect the hypervariable regions of the reference (_E. coli_) 16S. If the input file contains more
than one sequence it will output the prediction for all of them, and since this can be slow (and useless) it can be avoided with `-m INT` switch 
(maximum number of sequences to parse).


## Datasets
the `datasets` directory contains Illumina paired end sequences obtained from: whole genome sequencing, 16S amplicon sequencing, ITS amplicon sequencing.


