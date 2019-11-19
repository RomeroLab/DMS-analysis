#!/bin/bash

# set paths
data="/home/promero/illumina_data/PAR4/PAR4B/"
bt2="/home/promero/code/bowtie2-2.1.0/"


# prepare the six index files required by bowtie
"$bt2"bowtie2-build "$data"1GNXpet22_SgrAI_DraIII.fasta "$data"index_file
clear 


#BOWTIE2 options:
  #alignment types:
  # --end-to-end: requires ends of reads to match (default) - I like end-to-end because it is more stringent
  # --local: will trim off the ends of a read if it helps the alignment                                                                                                              


  # standard options:
  #--no-unal: don't write sequences that failed to align                                                       
  #--no-hd: don't have a SAM header                                                                                                                         
  #--maxins 2000: search for inserts that are up to 2000 bp (RE fragment was 1800 bp)
  # -p 7: use 7 cores
  #--dovetail: allow dovetails
  #--no-discordant: don't look for discordant alignments (instead just map reads individually)


  #--fr/--rf/--ff: relative orientation of reads (inward, outward, and same direction), default is inward


runs="Lib8v1_S1
Lib8v2_S2
Lib8v3_S3
Lib8pre_S4"

> alignment_output.txt

for r in $runs
do
    echo "$bt2"bowtie2 --very-sensitive-local --maxins 2000 --dovetail --no-discordant --no-unal --no-hd -p 7 -x "$data"index_file -1 "$data""$r"_L001_R1_001.fastq -2 "$data""$r"_L001_R2_001.fastq -S "$data""$r".sam 1>> alignment_output.txt
    time "$bt2"bowtie2 --very-sensitive-local --maxins 2000 --dovetail --no-discordant --no-unal --no-hd -p 7 -x "$data"index_file -1 "$data""$r"_L001_R1_001.fastq -2 "$data""$r"_L001_R2_001.fastq -S "$data""$r".sam 2>> alignment_output.txt
    echo >> alignment_output.txt
done

