import sequence_tools
import NGS_tools
from sys import argv


### input files and options ##########################################

samfile = argv[1] # the sam file to be converted to mutations
reffile = '1GNXpet22_SgrAI_DraIII.fasta' # the sequence used as a reference in the Bowtie2 run
outfile = samfile.split('/')[-1].replace('.sam','_mutations.txt')


start = 150 # the start and stop indices for the coding region.  Slicing reference[start:stop] should start on start codon and end on stop codon
stop = 1656

Qcut = 30 # Quality score cutoff (QS >= Qcut is considered good)

CDN_muts = False # return codon mutations (rather than AA mutations)
WT_reads = False # return reads that contain no mutations (necessary for counting WT to get frequencies) 



### load sequences and genetic code #################################

reference = sequence_tools.read_fasta(reffile)[1][0]
WTcodons = sequence_tools.split_codons(reference[start:stop])

code = sequence_tools.code
degcode = dict((p1+p2+'N',code[p1+p2+'A']) for p1 in 'ACGT' for p2 in 'ACGT' if len(set([code[p1+p2+p3] for p3 in 'ACGT']))==1) # these are the 8 AAs where the third position doen't matter
code.update(degcode)




### loop over all reads in sam file  #################################

output = ''
samfile = open(samfile)
sc = 0
lc = 0
while True:
    read = samfile.readline()
    if not read: break
    read_type = read.strip().split('YT:Z:')[1]

    if read_type=='CP': # concordant pair
        mate = samfile.readline()
        if not mate: break
        name1,align1 = NGS_tools.parse_SAM(read,reference)
        name2,align2 = NGS_tools.parse_SAM(mate,reference)
        if name1!=name2: print('Pair names do not match!',1/0)
        align1 = NGS_tools.remove_lowQ_indels(align1,Qcut)
        align2 = NGS_tools.remove_lowQ_indels(align2,Qcut)
        hiQ_indel = len([p for p in align1+align2 if '-' in p])>0 
        if hiQ_indel: continue # if there is a high quality indel, skip and go to the next read
        alignment = [align1[i] if align1[i][2]>align2[i][2] else align2[i] for i in range(len(align1))]

    elif read_type=='UP': #unpaired
        name,align = NGS_tools.parse_SAM(read,reference)
        alignment = NGS_tools.remove_lowQ_indels(align,Qcut)
        hiQ_indel = len([p for p in alignment if '-' in p])>0 
        if hiQ_indel: continue # if there is a high quality indel, skip and go to the next read


    # by this point there are no indels (first column of alignment==reference)
    alignment = alignment[start:stop]
    seqread = ''.join([p[1] if p[2]>=Qcut else 'N' for p in alignment])
    codons = sequence_tools.split_codons(seqread) 
    sc += 1 # add one to sequence count 

    # find the first and last codons that match the WT. Only use the codons between these. (this may help minimize reads into the adaptor sequences)
    matches = [i for i in range(len(codons)) if codons[i]==WTcodons[i]]
    if len(matches)==0: continue # nothing to do here


    # loop over all coding positions and find mutations
    mutations = []
    for i in range(len(codons)):
        if i>=matches[0] and i<=matches[-1]: # i is within matching range
            if codons[i] in code: # codes for a specific amino acid

                if CDN_muts: # store codon mutations
                    WT = WTcodons[i]
                    mut = codons[i]
                    if WT != mut:
                        mutations.append(WT+str(i)+mut)

                else: # store AA mutations 
                    WT = code[WTcodons[i]]
                    mut = code[codons[i]]
                    if WT != mut:
                        mutations.append(WT+str(i)+mut)

    if not WT_reads:
        if len(mutations)==0: continue

    mutations = ','.join(mutations)
    output += '%i; %i; %s\n' % (matches[0],matches[-1],mutations)
    lc += 1 # add one to line count 

open(outfile,'w').write(output)

print('Analyzed %i reads and wrote %i reads to file: %s' % (sc,lc,outfile))
