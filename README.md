# rnaseq
University Hospital Essen, NeuroScienceLab

In this repository I provide the code that I am using for RNA-sequencing data analysis in the course of my employment at the University Hospital Essen. We are working predominantly of bulk RNAseq data from either mice (middle cerebral artery occlusion model) or from human patients. The sequencing data is used for differential gene expression analysis and then for gene set enrichment analysis. 

I usually have long script from which I execute small code snippets consecutively, cause I find it convenient, fast and versatile. 

In my analysis I prepare the data with set of different tools for each step sometimes using more, sometimes fewer tools, eg:<br />
Trimming: raw / trimmomatic / fastp<br />
Mapping: star / hisat2 / tophat2 / salmon<br/>
Differential gene expression: DESeq2 / edgeR<br />
