#!/bin/bash
#SBATCH --job-name=crispr_count    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luyang@coh.org     # Where to send mail
#SBATCH -n 10                 # set the number of cores
#SBATCH -N -1
#SBATCH --mem=40G
#SBATCH --time=12:00:00               # Time limit hrs:min:sec, set max wallclock time
## SBATCH --output=crispr_count_%j.log   # Standard output and error log


####   Authored by Lu Yang, PhD, Department of System biology @ City of Hope.
#### v3 add sam to bed file conversion on 1/9/2018
#### v4 modify bowtie2 alignment options and downstream sam extract criterias on 1/10/2018
####   lib gRNA_id naming format:  Species_gene_pos_strand,  eg:  Ms_Dot1l_4129_-
#### v5a, add step to check duplicated sgRNA sequences and keep only one since it's all from the same gene with different synonyms,export duplicate list
#### v5b, add step to automatically check sgRNA length based on provided library, to the run_crispr_count_Apollo_v5.sh

start_time=$(date +%s.%N)       ## record script starting time
echo "Start crispr count running on hpc at $(date +'%r-%x')"

### start in root directory that contain fqfiles, crispr, bt2_index
flowID=$1             ## example: @NB501311
fqfile=$2             ## example: fqfiles/*.fq
pattern=$3            ## example: pattern = default or other gRNA structure pattern
bt2ref=$4             ## example: bt2_index/gRNA_Dot1l
gene=$5
nDup=$6
sgRNALen=$7

#echo "flowcell ID is $flowID"
#echo "current processed fastq file is $fqfile"
#echo "regular expression match pattern is $pattern"
printf "There are $nDup duplicated sgRNAs being removed in the provided library! ..................\n"
printf "Detected sgRNA length is $sgRNALen.\n"
echo "bowtie2 index path is $bt2ref"
rootDir=$PWD
echo "root directory is $rootDir"

# ## step 1: check lib files to see if there's any duplicated sequences, check sgRNA length
# cd bowtie_index
# refdir=$PWD
# LibName=$(basename $libfile)
# LibName=${LibName%.txt}
# #sgRNALen=$(head -n 1 $libfile |cut -f2|wc -c)
# sgRNALen=$(echo "$(head -n 1 $libfile |cut -f2|wc -m) - 1"|bc)
# printf "Detected sgRNA length is $sgRNALen...........\n"
# awk '!a[$2]++' $libfile > $LibName'_uniq.txt'
# n1=$(wc -l $libfile|cut -f1 -d" ")
# n2=$(wc -l $LibName'_uniq.txt'|cut -f1 -d" ")
# nUniq=$(echo "$n1 - $n2" | bc)
# printf "There are $nUniq duplicated sgRNAs being removed in the provided library! ..................\n"

# awk '{print ">"$1;print $2}' $LibName'_uniq.txt' > $LibName'.fa'
# bowtie2-build $LibName'.fa' $LibName >bt2_index.log 2>&1 
# bt2ref=$refdir/$LibName
# printf "Finished buidling bowtie2 references using unique sgRNAs library file at $(date +'%r')"
# cd ..

## step 1 extract gRNA sequences, working dir: fqfiles
cd fqfiles
echo "current working directory is $PWD"
fqdir=$PWD
fqfile=$(basename $fqfile)
#python $rootDir/scripts/gRNA_target.py $flowID $fqfile $pattern
#python /net/isi-dcnl/ifs/user_data/CWChen_Grp/Lu/testData_CRISPR/scripts/gRNA_target.py $flowID $fqfile $pattern
python /net/isi-dcnl/ifs/user_data/CWChen_Grp/Lu/pipelines/crispr/gRNA_target.py $flowID $fqfile $pattern

echo "finished gRNA target sequence extracting at $(date +'%r')"
cd ..

## step 2 map extracted reads to gRNA ref using bowtie2:
## 2a: bowtie2 -a option report all alignments
cd bowtie2_mapped
echo "current working directory is $PWD"
bt2dir=$PWD
#fname=${fqfile%.f*q}
fname=${fqfile%_L[0-9]*_R[0-9]_[0-9]*.fastq}
#bowtie2 -x $bt2ref -U $fqdir/$fname'_extracted.fastq' -S $fname'_extracted.sam' -p 4 --very-sensitive-local --time
## bowtie 2.3.4 if use very sensitive local option then no alignment will be reported, earlier version such as 2.1.0 is fine
bowtie2 -x $bt2ref -U $fqdir/$fname'_extracted.fastq' -S $fname'_extracted.sam' -a -p 4 --time
echo "finished mapping extracted gRNA target sequence using Bowtie2 at $(date +'%r')"
cd ..

#step 3: count reads perfect match to 20nt gRNA ref and summarized count table,
## exclude tag 16 = abadon reverse complement alignments, only keep tag == 0
## CIGAR 20M == perfect alignment match; MD tag: MD:Z:20 == perfect sequence match
cd crispr_count
echo "current working directory is $PWD"
echo -e "gRNA_id""\t"$fname > $fname'_extracted_count.txt'
perMpattern='MD:Z:'$sgRNALen
printf "$perMpattern \n"
awk '{if($2==0)print $0}' $bt2dir/$fname'_extracted.sam' |grep $perMpattern |awk '{print $3}' |sort|uniq -c |awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}'>> $fname'_extracted_count.txt'
#awk '{if($2==0)print $0}' $bt2dir/$fname'_extracted.sam' |grep MD:Z:20 |awk '{print $3}' |sort|uniq -c |awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}'>> $fname'_extracted_count.txt'
#awk '{if($2==0)print $0}' $bt2dir/$fname'_extracted.sam' |grep MD:Z:19 |awk '{print $3}' |sort|uniq -c |awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$1}'>> $fname'_extracted_count.txt'
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
rawreads=`echo $(wc -l $fqdir/$fqfile|cut -f1 -d" ")/4 | bc`   ## raw fastq reads number
rawreads1=$(echo $rawreads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta')    ## convert reads number to thousand separator
extractreads=`echo $(wc -l $fqdir/$fname'_extracted.fastq'|cut -f1 -d" ")/4 | bc`   ## raw fastq reads number
extractPer=$(echo `echo "scale=2; $extractreads*100/$rawreads" | bc`%)   ## gRNA structure extract percentage
samfile=$bt2dir/$fname'_extracted.sam'
totalalign=$(grep -v ^@ $samfile|wc -l)
PM=$(grep $perMpattern $samfile|wc -l)          # PM: Perfect match,  all 20nt match
PMper=$(echo `echo "scale=2; $PM*100/$totalalign" | bc`%)
# MM1=$(grep 19M $samfile|wc -l)         # MM1: one mismatch, 19M1S
# MM2=$(grep 18M $samfile|wc -l)         # MM2: two mismatch, 18M2S
# MM1per=$(echo `echo "scale=2; ($PM+$MM1)*100/$totalalign" | bc`%)
# MM2per=$(echo `echo "scale=2; ($PM+$MM1+$MM2)*100/$totalalign" | bc`%)
dur=$(echo "$(date +%s.%N) - $start_time" | bc)   ## calculate script execution time
dur1=$(echo "($(date +%s.%N) - $start_time)/60" | bc)
printf "Execution time: %.6f seconds\n" $dur
printf "Execution time: %.6f minutes\n" $dur1
## print to run_stats_table
#echo -e $fname"\t"$rawreads1"\t"$extractPer"\t"$PMper"\t"$MM1per"\t"$MM2per"\t"$dur >> run_stats.txt
echo -e $fname"\t"$rawreads1"\t"$extractPer"\t"$PMper"\t"$dur >> run_stats.txt


echo "Done crispr count running on hpc at $(date +'%r-%x')"
rm $fqdir/$fname'_extracted.fastq'