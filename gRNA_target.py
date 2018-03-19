### Function: extract 20nt gRNA target sequences from fastq files
### ARGS:     [keywords of flowcell ID][fastq files to open] [regular expression pattern to extract 20nt target sequence]
### ARGS examples:   @NB501311, B1.fq, default
### default pattern: CACCG[ACGT]{20}GTTTA
### Author:  Lu Yang @ system biology, city of hope, 12/25/2017


from sys import argv
script,flowID,fqfile,pattern =argv

## import required module
import re
import time
from datetime import timedelta     ### to benchmarking python programming time
start_time = time.time()


## verify ARGS
if pattern == "default":
    pattern = "CACCG([ACGT]{20})GTTTA"    ## the 20nt inside () is group 1 extracted info
    print "pattern used is %s" %(pattern)
else:
    print "pattern used is %s" %(pattern)
print "flowcell ID is "+flowID
print "fastq file ready to process is "+fqfile


###  open fastq, go through line by line, check if pattern exist in the current readid,
###  if match, extract the 20nt target sequence inside and export to new fastq output

# @NB501311:57:HYN75BGXY:1:11101:12416:1064 1:N:0:CCTCTGTA
# ATACANGATCTCTTGTGGAAAGGACGAAACTCCGCCTCTGGGATGGGAAGCTCAGTTTAAGAGCTATGCTGGAAA
# +
# AAAAA#EEEEEEEEE/EEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEAEAEEEEAAEEEE

#fname = re.sub(r".fastq|.fq","",fqfile)
fname = re.sub(r".fastq|.fq|_L[0-9]*_R1_[0-9]*.fastq","",fqfile)
with open(fqfile) as f:
    with open(fname+"_extracted.fastq", 'w') as fqoutput:
        for line in f:
            ### look for valid read id
            if line.startswith(flowID):
                readID = line.replace("\n","")
                readvalid = 1
                seqvalid = 0
                #print readID
            ### extract matched pattern target 20nt sequence and get index for the sequences
            elif re.search(pattern, line) and readvalid:
                match = re.search(pattern,line)
                sequence = match.group(1)
                seqstart = match.start(1)
                seqend = match.end(1)
                seqvalid = 1
            ### look for quality scores and use index from above to extract correspondent qual scores
            elif re.search(".{2,}",line) and seqvalid:
                qual = line[seqstart:seqend]
                print>>fqoutput,readID,'\n',sequence,'\n',"+",'\n',qual

print "Fastq data extracted successfully"

elapsed_time_secs = time.time() - start_time

print "Execution took: %s secs (Wall clock time)" % timedelta(seconds=round(elapsed_time_secs))
