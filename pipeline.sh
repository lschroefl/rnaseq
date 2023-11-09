#11/08/2023

#### ---------------- Architecure -----
# trimming: raw
# mapping: salmon / hisat2 (just for peace of mind regarding low mapping percentage) 
# differential gene expression: DESeq2
# gsea: clusterprofiler


# diverse renaming operations------------------------------
rename -n --verbose 's/Core/Penumbra/' *.fastq.gz
rename --verbose 's/Ipsilateral\-Penumbra\-trt/trt\-Ipsilaterl\-Penumbra/' *.fastq.gz
rename -n --verbose 's/^p/okcool/' *fastq.gz
rename --verbose 's/\_S[0-9]*/\-Ipsilateral\-Core/' c*.fastq.gz
rename --verbose 's/^c//' c*.fastq.gz
rename --verbose 's/_S[1-9]|10_/_/' 372*fastq.gz
rename --verbose 's/__/_/' 366*fastq.gz

# merging the files---------------------------------------------------
pathMain="/projects/neuronet/neuronet/processing"
pathOther="/projects/neuronet/neuronet/processing/raw/split"

sampleNr=($(seq 363 372))
sampleNr+=($(seq 416 430))
printf "%s\n" "${sampleNr[@]}"

cd $pathOther
for Nr in $(printf "%s\n" "${sampleNr[@]}")
do 
    if test -f $Nr-*-ipsiCore_L001_R1_001.fastq.gz; then
        file1=$(ls $Nr-*-ipsiCore_L001_R1_001.fastq.gz)
        file2=$(ls $Nr-*-ipsiPenumbra_L001_R1_001.fastq.gz)
        echo $file1
        cat $Nr-*-ipsiCore_L001_R1_001.fastq.gz >> ../concat/$file1
        cat $Nr-*-ipsiCore_L002_R1_001.fastq.gz >> ../concat/$file1
        cat $Nr-*-ipsiCore_L003_R1_001.fastq.gz >> ../concat/$file1
        cat $Nr-*-ipsiCore_L004_R1_001.fastq.gz >> ../concat/$file1
        
        echo $file2
        cat $Nr-*-ipsiPenumbra_L001_R1_001.fastq.gz >> ../concat/$file2
        cat $Nr-*-ipsiPenumbra_L002_R1_001.fastq.gz >> ../concat/$file2
        cat $Nr-*-ipsiPenumbra_L003_R1_001.fastq.gz >> ../concat/$file2
        cat $Nr-*-ipsiPenumbra_L004_R1_001.fastq.gz >> ../concat/$file2

    elif test -f $Nr-*-ipsiPenumbra_L001_R1_001.fastq.gz; then
        file3=$(ls $Nr-*-ipsiPenumbra_L001_R1_001.fastq.gz)
        echo $file3
        cat $Nr-*-ipsiPenumbra_L001_R1_001.fastq.gz >> ../concat/$file3
        cat $Nr-*-ipsiPenumbra_L002_R1_001.fastq.gz >> ../concat/$file3
        cat $Nr-*-ipsiPenumbra_L003_R1_001.fastq.gz >> ../concat/$file3
        cat $Nr-*-ipsiPenumbra_L004_R1_001.fastq.gz >> ../concat/$file3

    else
        echo "file number $Nr is unknown"
    fi
done
rename --verbose 's/_L001_R1_001//' *fastq.gz

#check differences of files--------------------------------
cmp split/427-day7-trt-ipsiPenumbra_L002_R1_001.fastq.gz concat/427-day7-trt-ipsiPenumbra.fastq.gz
diff -c split/416-day7-ctrl-ipsiCore_L001_R1_001.fastq.gz concat/416-day7-ctrl-ipsiCore_L001_R1_001.fastq.gz

# checking if number of reads is the same---------------------------------------
zcat ./concat/363-day3-ctrl-ipsiPenumbra.fastq.gz | grep @NB551534 | wc -l
zcat ./split/363-day3-ctrl-ipsiPenumbra_L00*_R1_001.fastq.gz | grep @NB551534 | wc -l

zcat ./concat/364-day3-ctrl-ipsiPenumbra.fastq.gz | grep @NB551534 | wc -l &&
zcat ./split/364-day3-ctrl-ipsiPenumbra_L00*_R1_001.fastq.gz | grep @NB551534 | wc -l

zcat ./concat/427-day7-trt-ipsiPenumbra.fastq.gz | grep @NB551534 | wc -l &&
zcat ./split/427-day7-trt-ipsiPenumbra_L00*_R1_001.fastq.gz | grep @NB551534 | wc -l

# quality control

ls raw/concat/*.fastq.gz | parallel -I % --jobs $(ls raw/concat/ | wc -l) fastqc --outdir qc %
multiqc --outdir qc qc/*

## quality control contamination 

ls /projects/neuronet/neuronet/processing/raw/concat/*.fastq.gz | parallel -I % --jobs $(ls raw/concat/*fastq.gz | wc -l) \
fastq_screen --conf /projects/neuronet/FastQ_Screen_Genomes/fastq_screen.conf \
--outdir /projects/neuronet/neuronet/processing/qc/contamination/ %

multiqc --outdir qc/contamination/ qc/contamination/*

## download reference genome and annotation, & checksum-----------------------------------
wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -O Mus_musculus.GRCm39.dna.primary_assembly.fa.gz 
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz -O Mus_musculus.GRCm39.110.gtf.gz

sum Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
sum Mus_musculus.GRCm39.110.gtf.gz

## building a transcriptome --------------------------------------------
gunzip -k Mus_musculus.GRCm39.dna.primary_assembly.fa.gz > Mus_musculus.GRCm39.dna.primary_assembly.fa &&
gunzip -k Mus_musculus.GRCm39.110.gtf.gz > Mus_musculus.GRCm39.110.gtf &&
gffread -w Mus_musculus.GRCm39.transcripts.primary_assembly.fa -g Mus_musculus.GRCm39.dna.primary_assembly.fa  Mus_musculus.GRCm39.110.gtf

gzip Mus_musculus.GRCm39.transcripts.primary_assembly.fa

rm Mus_musculus.GRCm39.dna.primary_assembly.fa &&
rm Mus_musculus.GRCm39.110.gtf &&
rm Mus_musculus.GRCm39.transcripts.primary_assembly.fa

# building a decoy aware index for salmon------------------------------------------
# kmer size changed to k < readlength/2

inputDirectory=genAnno/
outputDirectory=indices/salmon
inputFasta=Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
inputTranscripts=Mus_musculus.GRCm39.transcripts.primary_assembly.fa.gz

grep "^>" <(gunzip -c "$inputDirectory$inputFasta") | cut -d " " -f 1 > "$inputDirectory"Mus_musculus.GRCm39.salmon.decoys.txt
sed --in-place -e 's/>//g' "$inputDirectory"Mus_musculus.GRCm39.salmon.decoys.txt
cat "$inputDirectory$inputTranscripts" "$inputDirectory$inputFasta" > "$inputDirectory"Mus_musculus.GRCm39.salmon.gentrome.fa.gz

inputGentrome="$inputDirectory"Mus_musculus.GRCm39.salmon.gentrome.fa.gz
inputDecoys="$inputDirectory"Mus_musculus.GRCm39.salmon.decoys.txt

(salmon \
index -p 22 -k 23 \
-t "$inputGentrome" \
-d "$inputDecoys" \
-i "$outputDirectory" \
> "$outputDirectory"salmon.log 2>&1; \
echo -e "\nexit code $? " >> "$outputDirectory"salmon.log; \
timestamp=$(date +"%D %T"); echo "$timestamp" >> "$outputDirectory"salmon.log; \
) & disown -h %%

#### building an index for hisat2---------------------------------------
gunzip -k ./genAnno/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

(hisat2-build -p 24 ./genAnno/Mus_musculus.GRCm39.dna.primary_assembly.fa indices/hisat2/hisat2.index > indices/hisat2/hisat2.log 2>&1; 
echo exit code $? >> indices/hisat2/hisat2.log;
timestamp=$(date +"%D %T"); echo "$timestamp" >> indices/hisat2/hisat2.log;) &


# 9 QUANTIFYING WITH SALMON-------------------------------------------------
mapper=salmon
inputDir="raw/concat/"  
indexDir="indices/$mapper/"
outputFolder="mapping/$mapper/"

( for fastqgz in $(ls $inputDir*day14*.fastq.gz) 
do 
    outputPrefix="$(basename $fastqgz .fastq.gz)"
    outputDirectory="$outputFolder$outputPrefix"
    logFile="$outputDirectory.log"

    salmon quant --index "$indexDir" --libType A \
    -r $fastqgz -p 22 --validateMappings --seqBias \
    --gcBias -o $outputDirectory 2> "$logFile"; \
    echo exit code $? >> "$logFile"; 
done ) & disown -h %% 



### Mapping with hisat2----------------------------------------------------

mapper=hisat2
inputDirectory=./raw/concat/
( for fastqgz in $(ls $inputDirectory*day14*.fastq.gz) 
do 
    output="$(basename $fastqgz .fastq.gz)"
    output="./mapping/$mapper/$output"
    hisat2 --threads 15 --mm \
        -x ./indices/hisat2/hisat2.index \
        -U $fastqgz  \
        -S $output.sam 2> $output.log; \
        echo -e "\n\nexit code $?" >> $output.log
done ) & disown -h %%


# MAPPING STATISTICS------------------------------------------

fastqDirectory="./raw/concat"; 
inputDirectory="./mapping";
tsvFile="./mapping/statistics.tsv";
logFile="./mapping/statistisc.log";
tempDir="./mapping/tmp";

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
        id="$mapper/$sampleNumber"
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

# counting of the hisat2 aligned features------------------------------
# Without counting multimappers
gtfFile="./genAnno/Mus_musculus.GRCm39.110.gtf.gz"
mapper=hisat2
inputDirectory="./mapping/$mapper/"
logFile="$inputDirectory"featureCounts.log
(featureCounts \
    -g gene_id -T 20 -s 0 \
    --verbose \
    -a "$gtfFile" \
    -o "$inputDirectory"featureCounts \
    $(ls "$inputDirectory"*.sam) \
    > "$logFile" 2>&1; \
    echo -e "\nexit code $?" >> "$logFile";
    timestamp=$(date +"%D %T"); \
    echo -n "at the $timestamp" >> "$logFile";) \
    & disown -h %%

# transcript to gene map for salmon


gtfFile="./genAnno/Mus_musculus.GRCm39.110.gtf.gz"
tx2geneFile="./genAnno/Mus_musculus.GRCm39.110.tx2gene"


#ensembl gtf file, with gene names instead of ensembl id
zcat $gtfFile | awk -F "\t" '$3=="transcript"{split($9,a,";"); split(a[3],b,"\""); 
split(a[5],c,"\""); print b[2]","c[2]}' > "$tx2geneFile"


## CONTINUATION OF THE PIPLINE IN SCRIPT "flow_neuronet.Rmd"


rename --verbose 's/001/ipsiPenumbra/' *fastq.gz

