# rnaseq
University Hospital Essen, NeuroScienceLab

In this repository I provide the code that I am using for RNA-sequencing data analysis in the course of my employment at the University Hospital Essen. We are working predominantly of bulk RNAseq data from either mice (middle cerebral artery occlusion model) or from human patients. The sequencing data is used for differential gene expression analysis and then for gene set enrichment analysis. 

I usually have a long script from which I execute small code snippets consecutively, cause I find it convenient, fast and versatile. 

A small description to some of the scripts in this repository: 

**flow_neuronet.sh** -> a bash script that contains all the steps from receiving raw fastq files until the creation of a count table containing genes as rows and samples as columns(ie. quality control, mapping, some mapping statistics calculations).
**flow_neuronet.Rmd** -> a R markdown file that presents the results in html format, in a clean and neat way. My analysis continues with the count table (from flow_neuronet.sh), performs differential gene expression analysis, gene set enrichment anaylisis, clustering and diverse visualizations for an intuitive understanding and exploration of the data. 
transform_to_integer.pl -> a perl script that rounds the count table, in case there are floats. 

**alternative1.sh & alternative2.sh** -> a bash script that contains all the steps from receiving raw fastq files until the creation of a count table containing genes as rows and samples as columns, but this time comparing different 
- annotations (ucsc & ensembl),
- trimmers (fastp & trimmomatic),
- mappers (tophat2 & hisat2 & star & salmon & segemehl)
- DGE models (edgeR & DESeq2)
to see how much the choice of these tools influences the result of our differential gene expression analysis. Indeed the differences in the results between these tools is quite remarkable, and highlighted the variablility of resulst, the fundamental importance of tool choice, and further of the need of consistent tools and guidelines for RNAseq analysis.

**alternative1.R & alternative2.R** -> the R scripts belonging respectively to alternative1.sh & alternative2.sh2. They perform the differential expression analysis and the visualization of the results 
