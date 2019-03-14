# Bioinformatic Workshop: Analyzing ChIP-seq data

## ChIP-seq workshop guideline
- Part 1
    - Downloading a dataset from GEO.
    - FASTQ quality control.
    - Alignment.
    - Post-processing.
    - Peak calling.
    - Visualization with IGV.
- Part 2
    - Introduction to R and RStudio.
    - Dealing with `GenomicRanges`.
    - Differential analysis with DESeq2.

## Software
- [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- (?) [AfterQC](https://github.com/OpenGene/AfterQC)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Samtools](http://www.htslib.org/)
- [MACS2](https://github.com/taoliu/MACS)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [bedGraphToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig)
	- [fetchChromSizes](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes)
- [IGV](http://software.broadinstitute.org/software/igv/)
- [RStudio](https://www.rstudio.com/)
- [R](https://www.r-project.org/)
    - [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
    - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

## Data
GEO Dataset: [GSE112221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112221). PMID: [30136421](https://www.ncbi.nlm.nih.gov/pubmed/30136421)

__Samples SRAs__:

Sample names                | Experiment ID
----------------------------|------------------------
NL1_h3k27ac_h3k27ac_chipseq | SRR6880492
Cirr1_h3k27ac_chipseq       | SRR6880493    
HCC1_h3k27ac_chipseq        | SRR6880494    
Cirr3_h3k27ac_chipseq       | SRR6880495    
HCC3_h3k27ac_chipseq        | SRR6880496
NL1_input_chipseq           | SRR6880507, SRR6880508
Cirr1_input_chipseq         | SRR6880509, SRR6880510
HCC1_input_chipseq          | SRR6880511, SRR6880512
Cirr3_input_chipseq         | SRR6880513, SRR6880514
HCC3_input_chipseq          | SRR6880515, SRR6880516
