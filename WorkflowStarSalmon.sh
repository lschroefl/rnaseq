#checking the number of reads
for fastqgz in $(ls ./processing/raw/concatenated/)
do 
    echo concatenated/"$fastqgz" $(zcat processing/raw/concatenated/$fastqgz | echo $(($(wc -l)/4))) >> ./processing/raw/var/raw_read_counts.txt
done

# investigating why 4 files
zcat processing/raw/363-Ipsilateral-Core-ctrl_S1_L004_R1_001.fastq.gz | tail -n 4

#1 concatenate multiple fastq files
for i in {1..10}
do 
    sample_identifier=$(ls ./raw/unconcatenated/3*-Ipsilateral-Core-*_S"$i"_L004_R1_001.fastq.gz) 
    sample_identifier=${sample_identifier%-Ipsilateral-Core-*_S"$i"_L004_R1_001.fastq.gz*}
    sample_identifier=${sample_identifier##*"./raw/unconcatenated/"}
    echo $sample_identifier
    cat ./raw/unconcatenated/3*-Ipsilateral-Core-*_S"$i"_L001_R1_001.fastq.gz  >> ./raw/concatenated/"$sample_identifier".fastq.gz
    cat ./raw/unconcatenated/3*-Ipsilateral-Core-*_S"$i"_L002_R1_001.fastq.gz  >> ./raw/concatenated/"$sample_identifier".fastq.gz
    cat ./raw/unconcatenated/3*-Ipsilateral-Core-*_S"$i"_L003_R1_001.fastq.gz  >> ./raw/concatenated/"$sample_identifier".fastq.gz
    cat ./raw/unconcatenated/3*-Ipsilateral-Core-*_S"$i"_L004_R1_001.fastq.gz  >> ./raw/concatenated/"$sample_identifier".fastq.gz
done



# 2 download suited sham data set https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112348
fasterq-dump --outdir ./raw/concatenated/ SRR6899124 #or 
fastq-dump --outdir ./raw/concatenated/ --split-3 --gzip SRR6899123 & \
fastq-dump --outdir ./raw/concatenated/ --split-3 --gzip SRR6899124 & \
fastq-dump --outdir ./setraw/concatenated/ --split-3 --gzip SRR6899125 & 





# 3 quality control of my raw reads
>./raw/concatenated/qc/fastqc.log
ls raw/concatenated/*.fastq.gz | parallel -I % --jobs 12 fastqc --outdir raw/concatenated/qc/ % 2>> ./raw/concatenated/qc/fastqc.log
multiqc --outdir raw/concatenated/qc raw/concatenated/qc/* 2> ./raw/concatenated/qc/multiqc.log

# 4 clip the raw reads, with trimmomatic
ls raw/concatenated/*.fastq.gz | sort -n | parallel --jobs 12 -n 1 trimmomatic SE \
-summary clipped/trimmomatic/{/.}.summary -threads 5 {} clipped/trimmomatic/{/} \
SLIDINGWINDOW:4:20 MINLEN:25 TRAILING:3 2>> clipped/trimmomatic/trimmomatic.log

ls clipped/trimmomatic/*.fastq.gz | parallel -I % --jobs 12 fastqc --outdir ./clipped/trimmomatic/qc/ % 2>> ./clipped/trimmomatic/qc/fastqc.log
multiqc --outdir clipped/trimmomatic/qc/ clipped/trimmomatic/qc/* 2> ./clipped/trimmomatic/qc/multiqc.log

# 5 download reference genome + annotation

wget --timestamping \
'ftp://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz' \
-O annotation/mm39.fa.gz 
wget --timestamping \
'ftp://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz' \
-O annotation/mm39.refGene.gtf.gz


wget --timestamping \
'https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz' \
-O ensembl/GRCm39.dna.primary_assembly.fa.gz 
wget --timestamping \
'https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz' \
-O ensembl/GRCm39.109.gtf.gz

#checking for integrity via checksums
sum ensembl/GRCm39.dna.primary_assembly.fa.gz 
sum ensembl/GRCm39.109.gtf.gz


# 6 create bed file

gunzip annotation/mm39.refGene.gtf
bedparse gtf2bed annotation/mm39.refGene.gtf > annotation/mm39.refGene.bed

# 7 CREATE TRANSCRIPTOME FOR SALMON
gunzip -c genome/ensembl/GRCm39.109.gtf.gz > genome/ensembl/GRCm39.109.gtf & 
gunzip -c genome/ensembl/GRCm39.dna.primary_assembly.fa.gz > genome/ensembl/GRCm39.dna.primary_assembly.fa &
gffread -w ensembl/GRCm39.dna.primary_assembly.transcripts.fa -g ensembl/GRCm39.dna.primary_assembly.fa  ensembl/GRCm39.109.gtf
gzip ensembl/GRCm39.dna.primary_assembly.transcripts.fa


# 8.1 INDEX FOR STAR
# star does not accept gzipped files for index generation apparently (--readFilesCommand does not work)
genomeDir=indices/ensembl/star/
genomeFasta=./genome/ensembl/GRCm39.dna.primary_assembly.fa
genomeGTF=./genome/ensembl/GRCm39.109.gtf
outDir=./indices/ensembl/star/
logFile=./indices/ensembl/star/genomeGenerate.log

(STAR \
--runMode genomeGenerate \
--genomeDir "$genomeDir" \
--genomeFastaFiles "$genomeFasta" \
--sjdbGTFfile "$genomeGTF" \
--sjdbOverhang 54 \
--runThreadN 22 \
--outFileNamePrefix "$outDir" \
--genomeSAindexNbases 14 \
> "$logFile" 2>&1; \
echo -e "\n\nexit code $?" >> "$logFile"; \
timestamp=$(date +%T); echo " at $timestamp" >> "$logFile";\
) & disown -h %% 

#--sjdbOverhang adjusted to read length 
#--genomeSAindexNbases adjusted to genome size

# 8.2 INDEX FOR SALMON

# building a list of decoys and a file containing the transcriptome and the genome concatenated

inputDirectory=genome/ensembl/
outputDirectory=indices/ensembl/salmon/
inputFasta=GRCm39.dna.primary_assembly.fa
inputTranscripts=GRCm39.transcripts.primary_assembly.fa


grep "^>" <(gunzip -c "$inputDirectory$inputFasta") | cut -d " " -f 1 > "$inputDirectory"GRCm39.salmon.decoys.txt
sed --in-place -e 's/>//g' "$inputDirectory"GRCm39.salmon.decoys.txt
cat "$inputDirectory$inputTranscripts" "$inputDirectory$inputFasta" > "$inputDirectory"GRCm39.salmon.gentrome.fa.gz
 
inputGentrome="$inputDirectory"GRCm39.salmon.gentrome.fa
inputDecoys="$inputDirectory"GRCm39.salmon.decoys.txt

## salmon accepts gzipped files for index creation
# decoy aware index for salmon
# I change the kmer size to k < readlength/2

(
salmon \
index -p 22 -k 23 \
-t "$inputGentrome" -d \
"$inputDecoys" -i \
"$outputDirectory" \
> "$outputDirectory"salmon.log 2>&1; \
echo -e "\nexit code $? " >> "$outputDirectory"salmon.log; \
timestamp=$(date +"%D %T"); echo "$timestamp" >> "$outputDirectory"salmon.log; \
) & disown -h %%

# 9 QUANTIFYING WITH SALMON
genome=ensembl
clipper=trimmomatic
mapper=salmon
#inputDir="$clipper/concatenated/"
inputDir="clipped/$clipper/"  #alternative if clipper
indexDir="indices/$genome/$mapper/"
outputFolder="mapping/$genome/$clipper/$mapper/"

( for fastqgz in $(ls $inputDir*.fastq.gz) 
do 
    outputPrefix="$(basename $fastqgz .fastq.gz)"
    outputDirectory="$outputFolder$outputPrefix"
    logFile="$outputDirectory.log"

    salmon quant --index "$indexDir" --libType A \
    -r $fastqgz -p 33 --validateMappings --seqBias \
    --gcBias -o $outputDirectory 2> "$logFile"; \
    echo exit code $? >> "$logFile"; 
done ) & disown -h %% 


# 10 MAPPING WITH STAR

#before the run load the index into RAM (so it will NOT be loaded for each individual run): 

inputDirectory=raw/concatenated/
indexDir=indices/ensembl/star/
outputFolder=mapping/ensembl/raw/star/

genomeLoadDir=mapping/ensembl/raw/star/genomeLoad/
(
STAR \
--genomeDir "$indexDir" \
--outFileNamePrefix $genomeLoadDir \
--genomeLoad LoadAndExit \
> "$genomeLoadDir"LoadAndExit.log 2>&1 \
echo -e "\nexit code $?" >> "$genomeLoadDir"LoadAndExit.log; \
timestamp=$(date +%T); \
echo "at $timestamp" >> "$genomeLoadDir"LoadAndExit.log ; \
) &

# the actual mapping
( for fastqgz in $(ls $inputDirectory*.fastq.gz) 
do 
    outputPrefix="$(basename $fastqgz .fastq.gz)"
    outputDirectory="$outputFolder$outputPrefix/"
    logFile="$outputFolder$outputPrefix.log"
    mkdir $outputDirectory
    STAR \
    --genomeDir "$indexDir" \
    --runThreadN 22 \
    --readFilesIn $fastqgz \
    --outFileNamePrefix $outputDirectory \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM 30000000000 \
    --outSAMunmapped Within \
    > "$logFile" 2>&1; \
    echo -e "\nexit code $?" >> "$logFile"; \
    timestamp=$(date +%T); \
    echo -n " at $timestamp" >> "$logFile"; \
done ) & disown -h %%

# after the run remove the genome: 

STAR \
--genomeDir "$indexDir" \
--genomeLoad Remove \
--outFileNamePrefix $genomeLoadDir \
> "$genomeLoadDir"remove.log 2>&1 \
echo -e "\nexit code $?" >> "$genomeLoadDir"remove.log; \
timestamp=$(date +%T); \
echo "at $timestamp" >> "$genomeLoadDir"remove.log 

# 11 RENAMING STAR OUTPUT

mapper=star
for bam in $(ls mapping/ensembl/trimmomatic/$mapper/*/Aligned.sortedByCoord.out.bam)
do 
    newbam=${bam%*/Aligned.sortedByCoord.out.bam*}.bam
    mv $bam $newbam
done

# 12 MAPPING STATISTICS

genome=ensembl;
clipper=trimmomatic;
fastqDirectory="./clipped/$clipper"; 
inputDirectory="./mapping/$genome/$clipper";
tsvFile="./mapping/$genome/$clipper/statistics.tsv";
logFile="./mapping/$genome/$clipper/statistisc.log";
tempDir="./mapping/$genome/$clipper/tmp";

(  
mkdir -p "$tempDir";
echo -e "ID\treads\tsanity\tmapped\tmapped\tunmapped\tunmapped\tuniqmappers\tuniqmappers\tmultimappers\tmultimappers\tsplitmappers\tsplitmappers\tallmappings\tallmappings\tID" \
> "$tsvFile";
> "$logFile"; 
for mapper in $(ls "$inputDirectory")
do 
    for alignment in $(ls "$inputDirectory/$mapper/"*.[bs]am)
    do 
        timestamp=$(date +"%D %T")
        printf "$timestamp calculation for alignment "$alignment" started\n" >> "$logFile"
        #formulating the parameters for ID and for the sanity check
        sampleNumber=${alignment##*"$inputDirectory/$mapper/"};
        sampleNumber=${sampleNumber%.[bs]am*}; 
        fastqgz_mate1="$fastqDirectory/$sampleNumber"'.fastq.gz'; 
        echo "processing fastqgz_mate1: $fastqgz_mate1" >> "$logFile"
        id="$clipper/$mapper/$sampleNumber"
        read1=$(zcat "$fastqgz_mate1" | wc -l | xargs -n 1 bash -c 'echo $(($1 / 4))' args)
        mate1=$(samtools view -c -F 2304 $alignment)
        #calculating mapping statistics            
        echo $(samtools view -c -F 2304 $alignment | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > "$tempDir/reads.txt"  &
        echo $(samtools view -F 2308 "$alignment" |  wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > "$tempDir/mapped.txt" &
        printf "%.2f" $(samtools view  -F 2308 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > "$tempDir/mapped_ratio.txt" &
        echo $(samtools view -f 4 -F 2304 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > "$tempDir/unmapped.txt" &
        printf "%.2f" $(samtools view -f 4 -F 2304 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > "$tempDir/unmapped_ratio.txt" &
        echo $(samtools view -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > "$tempDir/uniqmappers.txt" &
        printf "%.2f" $(samtools view -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > "$tempDir/uniqmappers_ratio.txt" &
        echo $(samtools view -f 256  "$alignment" | cut -f 1 | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > "$tempDir/multimappers.txt" &
        printf "%.2f" $(samtools view -f 256 "$alignment" | cut -f 1 | sort | uniq | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > "$tempDir/multimappers_ratio.txt" &
        echo $(samtools view "$alignment" | awk 'match ($6, "N")' | cut -f 1 | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > "$tempDir/splitmappers.txt" &
        printf "%.2f" $(samtools view "$alignment" | awk 'match ($6, "N")' | cut -f 1 | sort | uniq | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > "$tempDir/splitmappers_ratio.txt" &
        echo $(samtools view -F 4 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > "$tempDir/mapped_total.txt" &
        printf "%.2f" $(samtools view -F 4 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > "$tempDir/mapped_total_ratio.txt" &
        wait;
        #adding the info from the temporary txt files to my table
        printf "$id\t\
        $(cat "$tempDir/reads.txt")\t\
        $(echo $(($mate1 == $read1)))\t\
        $(cat "$tempDir/mapped.txt")\t\
        $(cat "$tempDir/mapped_ratio.txt")\t\
        $(cat "$tempDir/unmapped.txt")\t\
        $(cat "$tempDir/unmapped_ratio.txt")\t\
        $(cat "$tempDir/uniqmappers.txt")\t\
        $(cat "$tempDir/uniqmappers_ratio.txt")\t\
        $(cat "$tempDir/multimappers.txt")\t\
        $(cat "$tempDir/multimappers_ratio.txt")\t\
        $(cat "$tempDir/splitmappers.txt")\t\
        $(cat "$tempDir/splitmappers_ratio.txt")\t\
        $(cat "$tempDir/mapped_total.txt")\t\
        $(cat "$tempDir/mapped_total_ratio.txt")\t\
        $id\n"\
        >> "$tsvFile";
        timestamp=$(date +"%D %T")
        printf "$timestamp calculation for alignment "$alignment" finished\n" >> "$logFile";
    done
done; 
) 2>> "$logFile" & disown -h %%


# 13 featurecounts for my star aligned data 
genome=ensembl
gtfFile="genome/$genome/GRCm39.109.gtf"
clipper=trimmomatic
mapper=star
inputDirectory="./mapping/$genome/$clipper/$mapper/"
logFile="$inputDirectory"featureCounts.log
(featureCounts \
    -g gene_id -T 20 -s 0 \
    --verbose \
    -a "$gtfFile" \
    -o "$inputDirectory"featureCounts \
    $(ls "$inputDirectory"*.bam) \
    > "$logFile" 2>&1; \
    echo -e "\nexit code $?" >> "$logFile";
    timestamp=$(date +"%D %T"); \
    echo -n "at the $timestamp" >> "$logFile";) \
    & disown -h %%

# 13.1 counting with htseq-count

genome=ensembl
gtfFile="genome/$genome/GRCm39.109.gtf"
clipper=raw
mapper=star
inputDirectory="./mapping/$genome/$clipper/$mapper/"
logFile="$inputDirectory"htseq.log

(
htseq-count \
-f bam \
--stranded=no \
--type=exon \
--idattr=gene_id \
--mode=union \
--with-header \
$(ls "$inputDirectory"*.bam) \
"$gtfFile" \
> "$inputDirectory"htseqCount \
2> "$logFile"; \
echo -e "\nexit code $?" >> "$logFile"; \
timestamp=$(date +"%D %T"); \
echo -n " datetime: $timestamp" >> "$logFile";
) & 



# 14 extract the crucial data from counts (featureCounts and htseq-counts)

inputDirectory=./mapping/ensembl/trimmomatic/star/;
output="$inputDirectory"featureCounts_countdata;
cut -f 1,7- "$inputDirectory"featureCounts | sed "s/\.\/mapping\/raw\/star\///g" | sed '1,2 s/.[bs]am//g' > $output


inputDirectory=mapping/ensembl/trimmomatic/star/;
input="$inputDirectory"htseqCount;
output="$inputDirectory"htseqCount_countdata;
sed 's/\.\/mapping\/ensembl\/trimmomatic\/star\///g' $input | sed '1,2 s/.[bs]am//g' > $output



# 15 rounding my counts from featurecounts (DESeq expects integer values)

# --> see script tranform_to_integers.pl

# 15 DESeq analysis of the star alignment see flow_may.R

# 16 extract the significant results from the table in R
awk -F "," 'NR>1 && $6<0.01 {print $1"\t"$3"\t"$6"\t"}' results/star_deseq2_ctr_vs_trt.csv > results/star_deseq2_ctr_vs_trt.sig
# no significant results.... I try without the sham animals

# 17 creating a gene-to-transcipt map for the tximport of the salmon mapped data

gtfFile2="./genome/ucsc/mm39.refGene.gtf"
gtfFile="./genome/ensembl/GRCm39.109.gtf"
tx2geneFile="./genome/ensembl/tx2gene_GRCm39.109.csv"
#ucsc gtf file, with gene names instead of gene id's
awk -F "\t" '$3=="transcript"{split($9,a,";"); split(a[2],b,"\""); 
split(a[1],c,"\""); print b[2]","c[2]}' "$gtfFile" > "$tx2geneFile"

#ensembl gtf file, with gene names instead of ensembl id
awk -F "\t" '$3=="transcript"{split($9,a,";"); split(a[3],b,"\""); 
split(a[5],c,"\""); print b[2]","c[2]}' "$gtfFile" > "$tx2geneFile"







## to do


# 1 look up fastq files to identify lane sample etc
zcat raw/unconcatenated/368-Ipsilateral-Core-trt_S6_L003_R1_001.fastq.gz | head
zcat raw/concatenated/368.fastq.gz | grep "@NB551534:173:HJ3FHAFX5:3:11401:25687:1048 1:N:0:6"
zcat clipped/trimmomatic/trimmomatic/369.fastq.gz | grep "@NB551534:173:HJ3FHAFX5:3:11401:25687:1048 1:N:0:6"



# 3 do featurecounts without fractions, without multimappers



# NOT URGENT ----------------------------------------------------


# check if i have metadata in my summarized experiment dataset
# isualize with igv -- POSTPONED, FIRST THE DEG analysis


#----------- - - - - -- 



