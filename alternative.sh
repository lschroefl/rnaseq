
# gzipping the fastq files
ls /projects/neuronet/sadler_liesz/fastq/*.fastq | parallel -I % --jobs 7 gzip %
##-> next time just download again that might be quicker

# counting the number of reads
zcat clipped/trimmomatic/SRR9130101_1.fastq.gz | echo $(($(wc -l)/4))
#checking for flow cell, lane, tile etc
for fastq in $(ls raw/*fastq.gz)
do
zcat $fastq | head -n 1
done
#check some stuff in annotation 
zcat annotation/mm10.fa.gz | grep -w "chr2" -B 50
#check for the numbers of features in the gtf file 
awk '{a[$3]++}END{for(k in a){print k,a[k]}}' annotation/mm10.refGene.gtf
#check number of transcripts in my new file
grep ">" -c genome/transcripts.mm10.fa 
## easy renaming of mulpt            echo $(samtools view -f 67 -F 256 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/primary.txt &
ile files
for file in $(ls genome/bowtie2/*id*)
do 
    mv $file ${file//.idx/}
done
# semicolon operator for sequentially execution, $? for exit code
#gunzipping some files
ls clipped/fastp/test/* | parallel -I % --jobs 12 gunzip %
# remove all the files with some pattern
find -name *.counts* | xargs rm


#fastq counts
{ zcat clipped/trimmomatic/SRR9130105_1.fastq.gz | wc -l | xargs -n 1 bash -c 'echo $(($1 / 4))' args; zcat clipped/trimmomatic/SRR9130105_2.fastq.gz | wc -l | xargs -n 1 bash -c 'echo $(($1 / 4))' args; } &
#read a specific line
samtools view mapping/fastp/hisat2/SRR9130101.sam | sed '1!d'
# how many perfect alignments 
samtools view mapping/fastp/hisat2/SRR9130101.sam | grep 'NM:i:0' | wc -l
#single mapped entries
samtools view mapping/fastp/hisat2/SRR9130101.sam | grep -Eo "NH:i:[0-9]+" | cut -d ":" -f 3 | sort -n | uniq -c  
#just a simple head 
samtools view mapping/fastp/hisat2/SRR9130101.sam | head
#simple line count divided by 2 
samtools view -F 256 mapping/fastp/hisat2/SRR9130101.sam | wc -l
# check how many reads in the sam file
samtools view -c -F 0x900 mapping/trimmomatic/star/SRR9130105.bam 






cd /projects/neuronet/sadler2/processing

#1 fastqc of raw data
ls raw/*.fastq.gz | parallel -I % --jobs 12 fastqc --outdir raw/qc %
#2 multiqc 
multiqc -outdir raw/qc raw/qc/*
#3 trimming with fastp 
ls raw/*.fastq.gz | sort -n | parallel --jobs 6 -n 2 \
fastp --in1 {1} --out1 clipped/fastp/{1/} \
--in2 {2} --out2 clipped/fastp/{2/} -3 20 -l 25 \
--json clipped/fastp/{1/}.report.json --html clipped/fastp/{1/}.report.html \
--thread 10
#4 trimming with trimmomatic
ls raw/*.fastq.gz | sort -n | parallel --jobs 6 -n 2 trimmomatic PE -summary \
clipped/trimmomatic/{1/}.summary -threads 9 {1} {2} \
clipped/trimmomatic/{1/} clipped/trimmomatic/{1/}.unpaired \
clipped/trimmomatic/{2/} clipped/trimmomatic/{2/}.unpaired \
ILLUMINACLIP:/homes/lschroefl/.conda/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
SLIDINGWINDOW:4:20 MINLEN:25 TRAILING:3
#5 Quality control for clipped reads
ls clipped/trimmomatic/*.fastq.gz | parallel --jobs 12 fastqc --threads 5 --outdir clipped/trimmomatic/qc
ls clipped/fastp/*.fastq.gz | parallel -I % --jobs 12 fastqc --threads 5 --outdir clipped/fastp/qc %
multiqc --outdir clipped/trimmomatic/qc clipped/trimmomatic/qc/*
multiqc --outdir clipped/fastp/qc clipped/fastp/qc/*
#6 downloading genome and annotation (mm10 & refSeq annotation as stated by sadler and liesz)
wget --timestamping \
'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz' \
-O annotation/mm10.fa.gz
wget --timestamping \
'ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz' \
-O annotation/mm10.refGene.gtf.gz
#7 create bed file -- bedops gtf2bed did NOT work for me -- bedparse gtf2bed did the job
gunzip annotation/mm10.refGene.gtf.gz
bedparse gtf2bed annotation/mm10.refGene.gtf > annotation/mm10.refGene.bed
#8 genome index for star
gunzip -c genome/mm10.fa.gz > genome/mm10.fa

STAR --runMode genomeGenerate --genomeDir genome/star/ \
--genomeFastaFiles genome/mm10.fa --runThreadN 60 \
--outFileNamePrefix genome/star/ --genomeSAindexNbases 14 

#9 creating transcript index for salmon
gffread -w genome/transcripts.mm10.fa -g genome/mm10.fa annotation/mm10.refGene.gtf
gzip genome/transcripts.mm10.fa

salmon index -k 31 -t genome/transcripts.mm10.fa.gz -i genome/salmon/ 2> genome/salmon/cmd.log

#10 creating index for segemehl
segemehl.x -x genome/segemehl/mm10.idx -d genome/mm10.fa 2> genome/segemehl/cmd.log

#11 creating index for tophat2
bowtie2-build genome/mm10.fa genome/bowtie2/mm10.idx > genome/bowtie2/cmd.log
bowtie2-build --large-index genome/mm10.fa genome/bowtie2/mm10.idx > genome/bowtie2/cmdLargeIndex.log

#12 creating index for hisat2
{ hisat2-build -p 22 genome/bowtie2/mm10.fa genome/hisat2/mm10 1> genome/hisat2/mm10.stdout.log 2> genome/hisat2/mm10.stderr.log; echo exit code $? >> genome/hisat2/mm10.stderr.log; } &

#13 mapping star ## NOT SUITED TO BE DONE VIA PARALLEL -- 30 Gb OF RAM PER JOB (loading of the genomes)
ls clipped/fastp/*.fastq.gz | sort -n | parallel --jobs 6 -n 2 STAR --genomeDir genome/star --runThreadN 10 \
--readFilesIn clipped/fastp/{1/} clipped/fastp/{2/} --outFileNamePrefix mapping/fastp/star/{2/.}. \
--outSAMattributes All --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat 

#13 mapping star, found a way to load the genome only once... not for each run
#before the run: 
mapper=star
clipper=raw
inputDirectory=clipped/$clipper/
for fastqgz in $(ls $inputDirectory*1.fastq.gz) 
do 
    mate1=$fastqgz
    mate2=${fastqgz%1*}2${fastqgz##*1}
    outputPrefix="$(basename $fastqgz _1.fastq.gz)"
    outputDirectory=mapping/$clipper/$mapper/
done

STAR --genomeDir genome/star \
--outFileNamePrefix $outputDirectory"genomeLoad." \
--genomeLoad LoadAndExit \
> $outputDirectory/genomeLoad.log 2>&1; \
echo -e "\nexit code $?" >> $outputDirectory/genomeLoad.log; \
timestamp=$(date +%T); \
echo "at $timestamp" >> $outputDirectory/genomeLoad.log;

#the STAR run

(for fastqgz in $(ls $inputDirectory*1.fastq.gz) 
do 
    mate1=$fastqgz
    mate2=${fastqgz%1*}2${fastqgz##*1}
    outputPrefix="$(basename $fastqgz _1.fastq.gz)"
    outputDirectory=mapping/$clipper/$mapper/
    STAR --genomeDir genome/$mapper --runThreadN 22 \
    --readFilesIn $mate1 $mate2 \
    --outFileNamePrefix $outputDirectory$outputPrefix \
    --outSAMattributes All --outSAMtype BAM \
    SortedByCoordinate --readFilesCommand zcat \
    --genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 \
    --outSAMunmapped Within KeepPairs \
    > mapping/$clipper/$mapper/$outputPrefix.log 2>&1;
    echo -e "\nexit code $?" >> mapping/$clipper/$mapper/$outputPrefix.log; \
    timestamp=$(date +%T); \
    echo "at $timestamp" >> mapping/$clipper/$mapper/$outputPrefix.log; \
done) & 

# when done: 
STAR --genomeDir genome/star \
--outFileNamePrefix $outputDirectory"genomeRemove." \
--genomeLoad Remove > mapping/raw/star/genomeRemove.log 2>&1


#14 mapping tophat --- SHOULD HAVE SAVED STDERROR IN LOG FILE
# bowtie wants to have the genome fa file in the same directory as the index
mapper=tophat2
clipper=trimmomatic
inputDirectory=clipped/$clipper/
{ for fastqgz in $(ls $inputDirectory*1.fastq.gz) 
do 
    mate1=$fastqgz
    mate2=${fastqgz%1*}2${fastqgz##*1}
    outputPrefix="$(basename $fastqgz _1.fastq.gz)"
    outputDirectory=mapping/$clipper/$mapper/$outputPrefix
    tophat2 -p 22 -o $outputDirectory genome/bowtie2/mm10 $mate1 $mate2 2> $outputDirectory.log; \
    echo exit code $? >> $outputDirectory.log;
done } & disown -h %%

#merging the bam output and the unmapped.bam file (for the maooing statistics)

samtools merge -u SRR9130102.bam SRR9130102/SRR9130102.bam SRR9130102/unmapped.bam &
samtools merge -u SRR9130103.bam SRR9130103/SRR9130103.bam SRR9130103/unmapped.bam &
samtools merge -u SRR9130104.bam SRR9130104/SRR9130104.bam SRR9130104/unmapped.bam &
samtools merge -u SRR9130105.bam SRR9130105/SRR9130105.bam SRR9130105/unmapped.bam &
samtools merge -u SRR9130106.bam SRR9130106/SRR9130106.bam SRR9130106/unmapped.bam &
samtools merge -u SRR9130101.bam SRR9130101/SRR9130101.bam SRR9130101/unmapped.bam &




#15 mapping with salmon
clipper=trimmomatic
inputDirectory=clipped/$clipper/
mapper=salmon
{ for fastqgz in $(ls $inputDirectory*1.fastq.gz) 
do 
    mate1=$fastqgz
    mate2=${fastqgz%1*}2${fastqgz##*1}
    outputPrefix="$(basename $fastqgz _1.fastq.gz)"
    outputDirectory=mapping/$clipper/$mapper/$outputPrefix
    salmon quant --index genome/salmon/ --libType A \
    -1 $mate1 -2 $mate2 -p 20 --validateMappings --seqBias \
    --gcBias -o $outputDirectory 2> mapping/$clipper/$mapper/$outputPrefix.cmd.log; \
    echo exit code $? >> mapping/$clipper/$mapper/$outputPrefix.cmd.log; 
done } &


#16 mapping with segemehl, fastp clipped fastq.gz had to be gunzipped for segemehl to understand
clipper=fastp
mapper=segemehl
inputDirectory=clipped/$clipper/
{ for fastqgz in $(ls $inputDirectory*1.fastq.gz) 
do 
    mate1=$fastqgz
    mate2=${fastqgz%1*}2${fastqgz##*1} 
    outputPrefix="$(basename $fastqgz _1.fastq.gz)"
    outputDirectory=mapping/$clipper/$mapper/$outputPrefix
    segemehl.x -t 33 -i genome/segemehl/mm10.idx \
    -d genome/bowtie2/mm10.fa \
    -q $mate1 -p $mate2 > $outputDirectory.sam 2> $outputDirectory.log; \
    echo -e "\n\nexit code $?" >> $outputDirectory.log;
done } & disown -h %%



#17 mapping with hisat2.. 
clipper=trimmomatic
mapper=hisat2
inputDirectory=clipped/$clipper/
{ for fastqgz in $(ls $inputDirectory*1.fastq.gz) 
do 
    mate1=$fastqgz
    mate2=${fastqgz%1*}2${fastqgz##*1}
    outputPrefix="$(basename $fastqgz _1.fastq.gz)"
    outputDirectory=mapping/$clipper/$mapper/$outputPrefix
    hisat2 --threads 22 --mm -x genome/hisat2/mm10 \
    -1 $mate1 -2 $mate2  \
    -S $outputDirectory.sam 2> $outputDirectory.log; \
    echo -e "\n\nexit code $?" >> $outputDirectory.log
done } &



#18 renaming files that are stupidly named 
clipper=raw
mapper=star
for bam in $(ls mapping/$clipper/$mapper/*Aligned.sortedByCoord.out.bam)
do 
    newbam=${bam%Aligned.sortedByCoord.out.bam*}.bam
    mv $bam $newbam
done

clipper=raw
mapper=tophat2
for bam in $(ls mapping/$clipper/$mapper/*/accepted_hits.bam)
do 
    newbam=${bam%/accepted_hits.bam*}.bam
    mv $bam $newbam
done

#19 BASIC MAPPING STATISTICS ---------- WORKING VERSION ## add ./clipped/raw folder for statistics
#would be nice to count alignments with more than 4 errors

#FRAGMENTS 
( 
echo -e "ID\tfragments\tsanity\tmapped\tmapped\tunmapped\tunmapped\tuniqmappers\tuniqmappers\tmultimappers\tmultimappers\tmappingstotal\tmappingstotal" > mapping/statistics/statistics_fragments.csv;
> mapping/statistics/statistics_fragments.log;
for clipper in $(ls clipped/)
do  
    for mapper in $(ls mapping/"$clipper"/)
    do  
        for alignment in $(ls mapping/"$clipper"/"$mapper"/*.[bs]am)
        do 
            timestamp=$(date +%T)
            printf "$timestamp calculation for alignment "$alignment" started\n" >> mapping/statistics/statistics_fragments.log
            #formulating the parameters and values for sanity check
            sampleNumber=${alignment##*"mapping/$clipper/$mapper/"};
            sampleNumber=${sampleNumber%.*am*}; 
            fastqgz_mate1=clipped/"$clipper"/"$sampleNumber"'_1.fastq.gz'; 
            id="$clipper/$mapper/$sampleNumber"
            mate1=$(zcat "$fastqgz_mate1" | wc -l | xargs -n 1 bash -c 'echo $(($1 / 4))' args)
            fragments=$(samtools view -c -f 64 -F 2304 $alignment)

            #calculating the statistics for my alignments 
            echo $(samtools view -c -f 64 $alignment | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/fragments.txt  &
            echo $(samtools view -f 66 -F 2304 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/mapped.txt &
            printf "%.2f" $(samtools view -f 66 -F 2304 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$fragments" | bc -l) > mapping/statistics/tmp/ratio_mapped.txt &
            echo $(samtools view -f 64 -F 2306 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/unmapped.txt & # do the alignment with star again, with unmapped reads find same for tophat2
            printf "%.2f" $(samtools view -f 64 -F 2306 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$fragments" | bc -l) > mapping/statistics/tmp/ratio_unmapped.txt & #aint gonna work if sam file does not contain unmapped reads
            echo $(samtools view -f 66 -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/uniqmappers.txt &
            printf "%.2f" $(samtools view -f 66 -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | xargs -n 1 echo "$(wc -l $1)/$fragments" | bc -l) > mapping/statistics/tmp/ratio_uniqmappers.txt &
            echo $(samtools view -f 322 "$alignment" | cut -f 1 | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/multimappers.txt &
            printf "%.2f" $(samtools view -f 322 "$alignment" | cut -f 1 | sort | uniq | xargs -n 1 echo "$(wc -l $1)/$fragments" | bc -l) > mapping/statistics/tmp/ratio_multimappers.txt &
            echo $(samtools view -f 66 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/mappingstotal.txt &
            printf "%.2f" $(samtools view -f 66 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$fragments" | bc -l) > mapping/statistics/tmp/ratio_mappingstotal.txt &
            wait;
            
            printf "$id\t\
            $(cat mapping/statistics/tmp/fragments.txt)\t\
            $(echo $(($mate1 == $fragments)))\t\
            $(cat mapping/statistics/tmp/mapped.txt)\t\
            $(cat mapping/statistics/tmp/ratio_mapped.txt)\t\
            $(cat mapping/statistics/tmp/unmapped.txt)\t\
            $(cat mapping/statistics/tmp/ratio_unmapped.txt)\t\
            $(cat mapping/statistics/tmp/uniqmappers.txt)\t\
            $(cat mapping/statistics/tmp/ratio_uniqmappers.txt)\t\
            $(cat mapping/statistics/tmp/multimappers.txt)\t\
            $(cat mapping/statistics/tmp/ratio_multimappers.txt)\t\
            $(cat mapping/statistics/tmp/mappingstotal.txt)\t\
            $(cat mapping/statistics/tmp/ratio_mappingstotal.txt)\n"\
            >> mapping/statistics/statistics_fragments.csv
            timestamp=$(date +%T)
            printf "$timestamp calculation for alignment "$alignment" finished\n" >> mapping/statistics/statistics_fragments.log
        done
    done
done; 
) 2>> mapping/statistics/statistics_fragments.log & disown -h %%


# READS
( 
echo -e "ID\tleft_reads\tsanity\tmapped\tmapped\tunmapped\tunmapped\tuniqmappers\tuniqmappers\tmultimappers\tmultimappers\tsplitmappers\tsplitmappers\tallmappings\tallmappings\tID\tright_reads\tsanity\tmapped\tmapped\tunmapped\tunmapped\tuniqmappers\tuniqmappers\tmultimappers\tmultimappers\tsplitmappers\tsplitmappers\tallmappings\tallmappings\tID" \
> mapping/statistics/statistics_reads.csv;
> mapping/statistics/statistics_reads.log;
for clipper in $(ls clipped/)
do 
    for mapper in $(ls mapping/"$clipper"/)
    do 
        for alignment in $(ls mapping/"$clipper"/"$mapper"/*.[bs]am)
        do 
            timestamp=$(date +%T)
            printf "$timestamp calculation for alignment "$alignment" started\n" >> mapping/statistics/statistics_reads.log
            #accessing the clipped fastq files
            sampleNumber=${alignment##*"mapping/$clipper/$mapper/"};
            sampleNumber=${sampleNumber%.*am*}; 
            fastqgz_mate1=clipped/"$clipper"/"$sampleNumber"'_1.fastq.gz'; 
            fastqgz_mate2=clipped/"$clipper"/"$sampleNumber"'_2.fastq.gz'; 

            #formulating the parameters for ID and for the sanity check
            id="$clipper/$mapper/$sampleNumber"
            read1=$(zcat "$fastqgz_mate1" | wc -l | xargs -n 1 bash -c 'echo $(($1 / 4))' args)
            read2=$(zcat "$fastqgz_mate2" | wc -l | xargs -n 1 bash -c 'echo $(($1 / 4))' args)
            reads_total=$(echo "$read1+$read2" | bc )
            mate1=$(samtools view -c -f 64 -F 2304 $alignment)
            mate2=$(samtools view -c -f 128 -F 2304 $alignment)
            
            #calculating mapping statistics            
            echo $(samtools view -c -f 64 -F 2304 $alignment | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/left_reads.txt  &
            echo $(samtools view -f 64 -F 2308 "$alignment" |  wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/left_mapped.txt &
            printf "%.2f" $(samtools view -f 64 -F 2308 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > mapping/statistics/tmp/left_mapped_ratio.txt &
            echo $(samtools view -f 68 -F 2304 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/left_unmapped.txt &
            printf "%.2f" $(samtools view -f 68 -F 2304 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > mapping/statistics/tmp/left_unmapped_ratio.txt &
            echo $(samtools view -f 64 -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/left_uniqmappers.txt &
            printf "%.2f" $(samtools view -f 64 -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > mapping/statistics/tmp/left_uniqmappers_ratio.txt &
            echo $(samtools view -f 320  "$alignment" | cut -f 1 | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/left_multimappers.txt &
            printf "%.2f" $(samtools view -f 320 "$alignment" | cut -f 1 | sort | uniq | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > mapping/statistics/tmp/left_multimappers_ratio.txt &
            echo $(samtools view -f 64 "$alignment" | awk 'match ($6, "N")' | cut -f 1 | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/left_splitmappers.txt &
            printf "%.2f" $(samtools view -f 64 "$alignment" | awk 'match ($6, "N")' | cut -f 1 | sort | uniq | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > mapping/statistics/tmp/left_splitmappers_ratio.txt &
            echo $(samtools view -f 64 -F 4 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/left_mapped_total.txt &
            printf "%.2f" $(samtools view -f 64 -F 4 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate1" | bc -l) > mapping/statistics/tmp/left_mapped_total_ratio.txt &
            echo $(samtools view -c -f 128 -F 2304 $alignment | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/right_reads.txt  &
            echo $(samtools view -f 128 -F 2308 "$alignment" |  wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/right_mapped.txt &
            printf "%.2f" $(samtools view -f 128 -F 2308 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate2" | bc -l) > mapping/statistics/tmp/right_mapped_ratio.txt &
            echo $(samtools view -f 132 -F 2304 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/right_unmapped.txt &
            printf "%.2f" $(samtools view -f 132 -F 2304 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate2" | bc -l) > mapping/statistics/tmp/right_unmapped_ratio.txt &
            echo $(samtools view -f 128 -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/right_uniqmappers.txt &
            printf "%.2f" $(samtools view -f 128 -F 2304 "$alignment" | grep --word-regexp "NH:i:1" | xargs -n 1 echo "$(wc -l $1)/$mate2" | bc -l) > mapping/statistics/tmp/right_uniqmappers_ratio.txt &
            echo $(samtools view -f 384 "$alignment" | cut -f 1 | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/right_multimappers.txt &
            printf "%.2f" $(samtools view -f 384 "$alignment" | cut -f 1 | sort | uniq | xargs -n 1 echo "$(wc -l $1)/$mate2" | bc -l) > mapping/statistics/tmp/right_multimappers_ratio.txt &
            echo $(samtools view -f 128 "$alignment" | awk 'match ($6, "N")' | cut -f 1 | sort | uniq | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/right_splitmappers.txt &
            printf "%.2f" $(samtools view -f 128 "$alignment" | awk 'match ($6, "N")' | cut -f 1 | sort | uniq | xargs -n 1 echo "$(wc -l $1)/$mate2" | bc -l) > mapping/statistics/tmp/right_splitmappers_ratio.txt &
            echo $(samtools view -f 128 -F 4 "$alignment" | wc -l | LC_ALL=en_US.UTF-8 awk '{ printf("%'"'"'d\n", $0) }') > mapping/statistics/tmp/right_mapped_total.txt &
            printf "%.2f" $(samtools view -f 128 -F 4 "$alignment" | xargs -n 1 echo "$(wc -l $1)/$mate2" | bc -l) > mapping/statistics/tmp/right_mapped_total_ratio.txt &
            wait;
    
            printf "$id\t\
            $(cat mapping/statistics/tmp/left_reads.txt)\t\
            $(echo $(($mate1 == $read1)))\t\
            $(cat mapping/statistics/tmp/left_mapped.txt)\t\
            $(cat mapping/statistics/tmp/left_mapped_ratio.txt)\t\
            $(cat mapping/statistics/tmp/left_unmapped.txt)\t\
            $(cat mapping/statistics/tmp/left_unmapped_ratio.txt)\t\
            $(cat mapping/statistics/tmp/left_uniqmappers.txt)\t\
            $(cat mapping/statistics/tmp/left_uniqmappers_ratio.txt)\t\
            $(cat mapping/statistics/tmp/left_multimappers.txt)\t\
            $(cat mapping/statistics/tmp/left_multimappers_ratio.txt)\t\
            $(cat mapping/statistics/tmp/left_splitmappers.txt)\t\
            $(cat mapping/statistics/tmp/left_splitmappers_ratio.txt)\t\
            $(cat mapping/statistics/tmp/left_mapped_total.txt)\t\
            $(cat mapping/statistics/tmp/left_mapped_total_ratio.txt)\t\
            $id\t\
            $(cat mapping/statistics/tmp/right_reads.txt)\t\
            $(echo $(($mate2 == $read2)))\t\
            $(cat mapping/statistics/tmp/right_mapped.txt)\t\
            $(cat mapping/statistics/tmp/right_mapped_ratio.txt)\t\
            $(cat mapping/statistics/tmp/right_unmapped.txt)\t\
            $(cat mapping/statistics/tmp/right_unmapped_ratio.txt)\t\
            $(cat mapping/statistics/tmp/right_uniqmappers.txt)\t\
            $(cat mapping/statistics/tmp/right_uniqmappers_ratio.txt)\t\
            $(cat mapping/statistics/tmp/right_multimappers.txt)\t\
            $(cat mapping/statistics/tmp/right_multimappers_ratio.txt)\t\
            $(cat mapping/statistics/tmp/right_splitmappers.txt)\t\
            $(cat mapping/statistics/tmp/right_splitmappers_ratio.txt)\t\
            $(cat mapping/statistics/tmp/right_mapped_total.txt)\t\
            $(cat mapping/statistics/tmp/right_mapped_total_ratio.txt)\n"\
            $id\t\
            >> mapping/statistics/statistics_reads.csv;
            timestamp=$(date +%T)
            printf "$timestamp calculation for alignment "$alignment" finished\n" >> mapping/statistics/statistics_reads.log;
        done
    done
done; 
) 2>> mapping/statistics/statistics_reads.log & disown -h %%


# 20 featureCounts to receive count matrix from my mappers

clipper=raw
for mapper in $(ls ./mapping/$clipper/)
do
    featureCounts -M --fraction -g gene_id -T 50 -p -s 0 --verbose\
        -a annotation/mm10.refGene.gtf \
        -o mapping/$clipper/$mapper/counts \
        mapping/$clipper/$mapper/SRR9130101.*am \
        mapping/$clipper/$mapper/SRR9130102.*am \
        mapping/$clipper/$mapper/SRR9130103.*am \
        mapping/$clipper/$mapper/SRR9130104.*am \
        mapping/$clipper/$mapper/SRR9130105.*am \
        mapping/$clipper/$mapper/SRR9130106.*am \
        2> mapping/$clipper/$mapper/featureCounts.log; \
        echo -e "\n\nexit code $?" >> mapping/$clipper/$mapper/featureCounts.log; \
        timestamp=$(date +%T); \
        echo " at $timestamp" >> mapping/$clipper/$mapper/featureCounts.log; \
done & disown -h %%



# 21 rename and transform my count table from featurecounts
for clipper in $(ls clipped/)
do 
    for mapper in $(ls mapping/"$clipper"/)
    do 
        cut -f 1,7- ./mapping/$clipper/$mapper/counts | sed "s/mapping\/$clipper\/$mapper\///g" | sed '1,2 s/.[bs]am//g' > ./mapping/$clipper/$mapper/countdataFloat
    done 
done 

clipper=raw
for mapper in $(ls mapping/"$clipper"/)
do 
    cut -f 1,7- ./mapping/$clipper/$mapper/counts | sed "s/mapping\/$clipper\/$mapper\///g" | sed '1,2 s/.[bs]am//g' > ./mapping/$clipper/$mapper/countdataFloat
done 



# 21 prepare salmon results to countmatrix, AINT NECESSARY i AM GONNA USE TXIMPORT

tar -czf mapping/fastp/salmon/all_samples.tar mapping/fastp/salmon/SRR9130101 mapping/fastp/salmon/SRR9130102 mapping/fastp/salmon/SRR9130103 mapping/fastp/salmon/SRR9130104 mapping/fastp/salmon/SRR9130105 mapping/fastp/salmon/SRR9130106

# 22 rounding all my countdata from featurecount with my perl script
perl transform_to_integer.pl

# 23 creating a transcript-to-gene map for the tximport of the salmon results

awk -F "\t" '$3=="transcript"{split($9,a,";"); split(a[2],b,"\""); 
split(a[1],c,"\""); print b[2]","c[2]}' annotation/mm10.refGene.gtf > annotation/tx2gene_mm10.csv




#TO DO sadler

# do the counting again 
# RUN STEP


# 1 check whats wrong with the basic mapping stats 
# star does only count 10 M input reads instead of 26 M... 
# star counts the input fragments not the reads in the log.final.out
# do star and tophat NOT list the reads that are unmapped? 



# 2 visualize alignments with igv
# 3 salmon tximeta import
# 4 DESeq


