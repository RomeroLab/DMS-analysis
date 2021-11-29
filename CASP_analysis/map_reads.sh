#!/bin/bash

# set paths
bt2="/usr/local/bin/" 
# find .fasta reference file, and generate list of fastq files
not_contained_in () { 
    local array="$1[@]"
    local seeking=$2
    local in=0
    for element in "${!array}"; do
        if [[ $element == "$seeking" ]]; then
            in=1
            break
        fi
    done
    return $in
}

data=$(pwd)
runs=( )

for file in "$data/fastq/"/*
do
	a="$file"
	xpath=${a%/*} 
	xbase=${a##*/}
	xfext=${xbase##*.}
	xpref=${xbase%.*}
	
	if [ "${xfext}" == "fastq" ] && not_contained_in runs "${xpref%????????????}"; then
    runs+=(${xpref%????????????})
  fi

  if [ "${xfext}" == "fasta" ]; then
    reference_fasta=${a}        
  fi
        	
	#echo;echo path=${xpath};echo pref=${xpref};echo ext=${xfext}	
done

echo
echo "These are the fastq files you'll be aligning:"
echo ${runs[*]}
echo
echo "Your reference fasta is:"
echo $reference_fasta

read -r -p "Gucci? [y/n] " response
response=${response,,}    # tolower
if [[ "$response" =~ ^(yes|y|yerp|yeah)$ ]]; then
# prepare the six index files required by bowtie

  "$bt2"bowtie2-build $reference_fasta "$data/fastq/"index_file
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

  > alignment_output.txt # make empty log file

  for r in ${runs[@]} #align all fastq files in runs to reference fasta
  do
        echo "$bt2"bowtie2 --very-sensitive-local --maxins 2000 --dovetail --no-discordant --no-unal --no-hd -p 16 -x "$data/fastq/"index_file -1 "$data/fastq/""$r"_L001_R1_001.fastq -2 "$data/fastq/""$r"_L001_R2_001.fastq -S "$data/fastq/""$r".sam 1>> alignment_output.txt
        time "$bt2"bowtie2 --very-sensitive-local --maxins 2000 --dovetail --no-discordant --no-unal --no-hd -p 16 -x "$data/fastq/"index_file -1 "$data/fastq/""$r"_L001_R1_001.fastq -2 "$data/fastq/""$r"_L001_R2_001.fastq -S "$data/fastq/""$r".sam 2>> alignment_output.txt
        echo >> alignment_output.txt
  done
  ### You may need to edit the file suffix that the above command searches for. 
  # move the newly created samfiles and log file to a new directory
  mkdir sam_files
  mv *.sam ../sam_files
  mv alignment_output.txt sam_files
  cp $reference_fasta sam_files

else
  echo 'fuckin figure it out: Make sure your fastq files are unique and you only have one reference fasta in this directory.'
fi

