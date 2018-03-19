#!/bin/bash
### working dir is root folder
### command example: bash scripts/run_crispr_count_hpc.sh @NB501311 default $PWD/bowtie_index/gRNA_Dot1l_v2 Dot1l

## load all required modules#####
module load Bowtie2            ##
#################################

flowID=$1
pattern=$2
libfile=$3
gene=$4
#scripts=/net/isi-dcnl/ifs/user_data/CWChen_Grp/Lu/testData_CRISPR/scripts
scripts=/net/isi-dcnl/ifs/user_data/CWChen_Grp/Lu/pipelines/crispr

if [ ! -d "Ref" ]; then
  echo "reference directory does not exist, please check library file"
  exit 1
else
  cd Ref
  refdir=$PWD
  LibName=$(basename $libfile)
  LibName=${LibName%.txt}
  #sgRNALen=$(head -n 1 $libfile |cut -f2|wc -c)
  sgRNALen=$(echo "$(head -n 1 $libfile |cut -f2|wc -m) - 1"|bc)
  printf "Detected sgRNA length is $sgRNALen...........\n"
  awk '!a[$2]++' $libfile > $LibName'_uniq.txt'
  awk 'a[$2]++ == 1 { print $1,$2 " is duplicated"}'  $libfile > $LibName'_dup.txt'   ## export duplication list
  cut -f2 -d" " $LibName'_dup.txt' > temp
  grep -f temp $libfile |sort -k2,2 >> $LibName'_dup.txt'
  rm temp
  n1=$(wc -l $libfile|cut -f1 -d" ")
  n2=$(wc -l $LibName'_uniq.txt'|cut -f1 -d" ")
  nDup=$(echo "$n1 - $n2" | bc)
  printf "There are $nDup duplicated sgRNAs being removed in the provided library! ..................\n"
  ## generate library files with unique sgRNA sequences
  awk '{print ">"$1;print $2}' $LibName'_uniq.txt' > $LibName'.fa'
  bowtie2-build $LibName'.fa' $LibName >bt2_index.log 2>&1 
  bt2ref=$refdir/$LibName
  printf "Finished buidling bowtie2 references using unique sgRNAs library file at $(date +'%r')\n"
  cd ..
fi

if [ ! -d "Apollo_logs" ]; then
  mkdir Apollo_logs
fi
if [ ! -d "crispr_count" ]; then
  mkdir crispr_count
else
  echo "crispr_count directory already exists, please remove then rerun the script"
  exit 1
fi
if [ ! -d "bowtie2_mapped" ]; then
  mkdir bowtie2_mapped
else
  echo "bowtie2_mapped directory already exists, please remove then rerun the script"
  exit 1
fi
if [ ! -d "bedfiles" ]; then
  mkdir bedfiles
else
  echo "bedfile directory already exists, please remove then rerun the script"
  exit 1
fi



batch=$(date +%Y%m%d%H%M%S)  ## record batch info by starting time eg 20180110093544
echo $batch
#touch run_stats.txt
if [ ! -f "run_stats.txt" ]; then
  #echo -e sampleID"\t"rawReads"\t""Extraction%""\t""PM%""\t""oneMM%""\t""twoMM%""\t""runtime(sec)" > run_stats.txt
  echo -e sampleID"\t"rawReads"\t""Extraction%""\t""PM%""\t""runtime(sec)" > run_stats.txt
else
  echo "run_stats.txt file exists"
  exit 1
fi






#$fname"\t"$rawreads"\t"$extractPer"\t"$PMper"\t"$MM1per"\t"$MM2per"\t"$dur
for i in fqfiles/*.f*q
do fname=${i%.fastq}
    fname=$(basename $fname)
	  echo $fname
		sbatch --output Apollo_logs/$fname'.crispr.count.'$batch'.log' $scripts/Slurm_crispr_count_v5.sh $flowID $i $pattern $bt2ref $gene $nDup $sgRNALen
done
