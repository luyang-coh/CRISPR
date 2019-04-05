# CRISPR

High density count pipeline
#Authored by Lu Yang, PhD, Department of System biology @ City of Hope.  
#v3 add sam to bed file conversion on 1/9/2018  
#v4 modify bowtie2 alignment options and downstream sam extract criterias on 1/10/2018  
#lib gRNA_id naming format: Species_gene_pos_strand, eg: Ms_Dot1l_4129_-    
#v5a, add step to check duplicated sgRNA sequences and keep only one since it's all from the same gene with different synonyms,export duplicate list   
#v5b, add step to automatically check sgRNA length based on provided library, to the run_crispr_count_Apollo_v5.sh   

#step 0: make ref fa file from gRNA_Dot1l_sequences check if duplicated sequences exist in    

awk '{print ">"$1;print $2}' $ref'.txt' > $ref'.fa' 
module load Bowtie2 
bowtie2-build $FA $NAME

#step 1: extract 20nt matched gRNA reads from original fastq files python gRNA_target.py @NB501311 *.fq pattern/default 
pattern: CACCG([ACGT]{20})GTTTA   (default);                                                                                                                               
User can also provided cutomized patterns for extractions


#step 2: map extracted reads to gRNA ref by bowtie2.

bowtie2 -x bt2.ref -U *_extracted.fastq -S *.sam -p 4 --very-sensitive-local --time


#step 3: count reads perfect match to 20nt gRNA ref and summarized count table 

grep 20M B7_C2_0_extracted.sam|awk '{print $3}' |sort|uniq -c |awk '{print $1,"\t",$2}'> *_count.txt

# Running on HPC suncluster
#current work dir: crispr root folder
#root folder structure: contains 
fqfiles ---- where raw fastq files store  
bowtie_index ---- gRNA lib bowtie index files 

EXAMPLE: 

bash /Your-Path-To/run_crispr_count_hpc.sh @NB501311 default $PWD/bowtie_index/sgKat8 Kat8

# Running on slurm server Apollo
bash /Your-Path-To/run_crispr_count_Apollo.sh @NB501311 default $PWD/bowtie_index/sgKat8 Kat8

QC pipeline
