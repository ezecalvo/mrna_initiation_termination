# mRNA initiation and termination is spatially coordinated
Ezequiel Calvo-Roitberg, Christine L. Carroll, GyeungYun Kim, Valeria Sanabria, Sergey V. Venev, Steven T. Mick, Joseph D. Paquette, Maritere Uriostegui-Arcos, Job Dekker, Ana Fiszbein, Athma A. Pai


### Short-read sequencing analysis

Short read sequencing data was processed and analyzed using the [HITindex pipeline](https://github.com/thepailab/HITindex) as outlined in the methods. Briefly, fastq files are mapped to genome assembly GRCh38.p14 for human or GRCm38.95 for mouse. The resulting bam files were used as input for the HITindex pipeline. First, an annotated metaexon file is generated from a gtf file (corresponding to the relevant genome assembly used to create the bam files) by collapsing overlapping consituent exons. Next, the pipeline uses the bam file and annotated metaexon file to to classify metaexons into an exon-type (.exon output) as well as calculating relative usage (PSI) of first and last exons (.AFEPSI and .ALEPSI outputs).The resulting AFEPSI and ALEPSI outputs from the HITindex pipeline are used for the short-read data analysis. While each RNA sequencing sample is handled individually within the HITindex, we load relevant AFEPSI and ALEPSI outputs into R as lists in order to perform bulk analysis across sample groups

### Long-read sequencing PITA genes identification
The _long_read_filtering_and_correlations.R_ script analyzes and filters long-read sequencing bed files after using [bam2bed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) split files where each feature of each read is one individual row. It outputs Spearman correlations between TSS-TES, filtered reads matching FE and polyA databases or any other databases that want to be checked.

Fastq files from PacBio reads were downloaded from the [ENCODE Project data portal](https://www.encodeproject.org/) and mapped to the hg38 human or mm38 reference genome using [minimap2](https://github.com/lh3/minimap2) following the developers’ recommended parameters: -ax splice:hq -uf and -ax splice. Reads were divided into multiple split features to define exons using [bedtools bam2bed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) and assigned to overlapping genes. Only primary alignments and reads with all their features assigned to the same gene were used as input for _long_read_filtering_and_correlations.R_.


### DRB4sU-LRSseq
The pythonscript_PITA_drb4suseq.py script identifies reads overlapping 5' end with annotated AFEs and 3’ end that do not overlap an annotated PAS peaks. These represent molecules that are still being synthesized by the RNA polymerase. For cDNA, it can also filter for presence of adapter in the 5' end to identify non-truncated reads, or on both 5' and 3' end adapters to identify complete reads.

The code requires as input a fastq.gz file containing all reads, coordinates of AFE (for example from first_exons_hit_index.bed), PAS peaks (for example from polyA_db.bed), and gene coordinates. For cDNA reads, it also requires a porechop output file.  The output is a info.gz file containing only the read features and if specified with the –outbam option, a .bam with the filtered reads. The downstream analysis was done with the LRS4sUDRB_Analysis.R script.
