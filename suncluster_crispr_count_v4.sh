#!/bin/bash
#$ -M luyang@coh.org
#$ -m bea
#$ -N crispr_count
#$ -pe shared 8
#$ -j yes
#$ -cwd
#$ -V

####   Authored by Lu Yang, PhD, Department of System biology @ City of Hope.
#### v3 add sam to bed file conversion on 1/9/2018
#### v4 modify bowtie2 alignment options and downstream sam extract criterias on 1/10/2018
####   lib gRNA_id naming format:  Species_gene_pos_strand,  eg:  Ms_Dot1l_4129_-,  sense: + ; anti-sense: -

start_time=$(date +%s.%N)       ## record script starting time
echo "Start crispr count running on hpc at $(date +'%r-%x')"

### start in root directory that contain fqfiles, crispr, bt2_index
flowID=$1             ## example: @NB501311
fqfile=$2             ## example: fqfiles/*.fq
pattern=$3            ## example: pattern = default or other gRNA structure pattern
bt2ref=$4             ## example: bt2_index/gRNA_Dot1l
gene=$5
echo "flowcell ID is $flowID"
echo "current processed fastq file is $fqfile"
echo "regular expression match pattern is $pattern"
echo "bowtie2 index path is $bt2ref"

###  could add command here to check if folder "crispr_count, bowtie2_mapped" exist already later
rootDir=$PWD
echo "root directory is $rootDir"

## step 1 extract gRNA sequences, working dir: fqfiles
cd fqfiles
echo "current working directory is $PWD"
fqdir=$PWD
fqfile=$(basename $fqfile)
#python $rootDir/scripts/gRNA_target.py $flowID $fqfile $pattern
python /isi-dcnl/ifs/user_data/CWChen_Grp/Lu/testData_CRISPR/scripts/gRNA_target.py $flowID $fqfile $pattern
echo "finished gRNA target sequence extracting at $(date +'%r')"
cd ..

## step 2 map extracted reads to gRNA ref using bowtie2:
## 2a: bowtie2 -a option report all alignments
## 2b: exclude tag 16 = abadon reverse complement alignments
cd bowtie2_mapped
echo "current working directory is $PWD"
bt2dir=$PWD
fname=${fqfile%.f*q}
#bowtie2 -x $bt2ref -U $fqdir/$fname'_extracted.fastq' -S $fname'_extracted.sam' -p 4 --very-sensitive-local --time
bowtie2 -x $bt2ref -U $fqdir/$fname'_extracted.fastq' -S $fname'_extracted.sam' -a -p 4 --very-sensitive-local --time
echo "finished mapping extracted gRNA target sequence using Bowtie2 at $(date +'%r')"
cd ..

#step 3: count reads perfect match to 20nt gRNA ref and summarized count table,
## exclude tag 16 = abadon reverse complement alignments, only keep tag == 0
## CIGAR 20M == perfect alignment match; MD tag: MD:Z:20 == perfect sequence match
cd crispr_count
echo "current working directory is $PWD"
echo -e "gRNA_id""\t"$fname > $fname'_extracted_count.txt'
#awk '{if($2==0)print $0}' $bt2dir/$fname'_extracted.sam' |grep MD:Z:20 |awk '{print $3}' |sort|uniq -c |awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}'>> $fname'_extracted_count.txt'
awk '{if($2==0)print $0}' $bt2dir/$fname'_extracted.sam' |grep MD:Z:19 |awk '{print $3}' |sort|uniq -c |awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}'>> $fname'_extracted_count.txt'
#grep MD:Z:20 $bt2dir/$fname'_extracted.sam'|awk '{print $3}' |sort|uniq -c |awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}'>> $fname'_extracted_count.txt'
echo "finished gRNA target sequence counting at $(date +'%r')"
cd ..

#step 4: generate bed files from count table for downstream peak calling: use cut position as start and extend only 1bp
# bed file format(0-based off): 1-chrname; 2-start; 3-end; 4-name; 5-score; 6-strand; (6 cols required by macs2)
# example of samfile column 3: Ms_Dot1l_2638_-
cd bedfiles
echo "current working directory is $PWD"
awk '{if($2==0)print $0}' $bt2dir/$fname'_extracted.sam' |grep MD:Z:20 |grep -E $gene |cut -f3|awk 'BEGIN{FS="_";OFS="\t"}{print $1 "_" $2,$3-1,$3,$2,0,$4}'> $fname'_extracted.bed'
echo "finished sam to bed file conversion at $(date +'%r')"
cd ..

##step 5: summarize run stats
rawreads=`echo $(wc -l $fqdir/$fname.fastq|cut -f1 -d" ")/4 | bc`   ## raw fastq reads number
rawreads1=$(echo $rawreads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta')    ## convert reads number to thousand separator
extractreads=`echo $(wc -l $fqdir/$fname'_extracted.fastq'|cut -f1 -d" ")/4 | bc`   ## raw fastq reads number
extractPer=$(echo `echo "scale=2; $extractreads*100/$rawreads" | bc`%)   ## gRNA strutcure extract percentage
samfile=$bt2dir/$fname'_extracted.sam'
totalalign=$(grep -v ^@ $samfile|wc -l)
PM=$(grep MD:Z:20 $samfile|wc -l)          # PM: Perfect match,  all 20nt match
MM1=$(grep 19M $samfile|wc -l)         # MM1: one mismatch, 19M1S
MM2=$(grep 18M $samfile|wc -l)         # MM2: two mismatch, 18M2S
PMper=$(echo `echo "scale=2; $PM*100/$totalalign" | bc`%)
MM1per=$(echo `echo "scale=2; ($PM+$MM1)*100/$totalalign" | bc`%)
MM2per=$(echo `echo "scale=2; ($PM+$MM1+$MM2)*100/$totalalign" | bc`%)
dur=$(echo "$(date +%s.%N) - $start_time" | bc)   ## calculate script execuation time
printf "Execution time: %.6f seconds" $dur
## print to run_stats_table
echo -e $fname"\t"$rawreads1"\t"$extractPer"\t"$PMper"\t"$MM1per"\t"$MM2per"\t"$dur >> run_stats.txt


echo "Done crispr count running on hpc at $(date +'%r-%x')"
