import sequence_tools
import pickle
import NGS_tools
from sys import argv
from time import time

# load the reference sequences used in the Bowtie2 run
reference = sequence_tools.read_fasta('caspase3_ref.fasta')[1][0]
coding = range(149,926)
WTcodons = sequence_tools.split_codons(''.join([reference[p] for p in coding]))
seqitr = xrange(len(WTcodons))
#1/0
# Define the quality cutoff for a good read
Qcut = 30 # 30 and above is good


# load amino acids and the genetic code
AAs = tuple(set(sequence_tools.code.values()))
code = sequence_tools.code
degcode = dict((p1+p2+'N',code[p1+p2+'A']) for p1 in 'ACGT' for p2 in 'ACGT' if len(set([code[p1+p2+p3] for p3 in 'ACGT']))==1) # these are the 8 AAs where the third position doen't matter
code.update(degcode)


# the dicts to keep counts
CDN_count = dict((i,dict((a,0) for a in code)) for i in range(len(coding)/3)) # dict of dicts
AA_count = dict((i,dict((a,0) for a in AAs)) for i in range(len(coding)/3)) # dict of dicts
pair_count = dict(((i,j),dict()) for i in range(len(coding)/3) for j in range(len(coding)/3) if i<j)


filename = argv[1]
samfile = open(filename)
PE_size = []
while True:
    read = samfile.readline()
    if not read: break
    read_type = read.strip().split('YT:Z:')[1] # CP for concordant pair, UP for unpaired
    if read_type=='CP':
        mate = samfile.readline()
        name1,align1 = NGS_tools.parse_SAM(read,reference)
        name2,align2 = NGS_tools.parse_SAM(mate,reference)
        if name1!=name2: print ('fuckup: pair names do not match!'), 1/0
        align1 = NGS_tools.remove_lowQ_indels(align1,Qcut)
        align2 = NGS_tools.remove_lowQ_indels(align2,Qcut)
        hiQ_indel = len([p for p in align1+align2 if '-' in p])>0 
        if hiQ_indel: continue # if there is a high quality indel, skip and go to the next read
        alignment = [align1[i] if align1[i][2]>align2[i][2] else align2[i] for i in range(len(align1))]

        # store some info about the size of the PE fragment
        obs_pos = [i for i,p in enumerate(alignment) if p[1]!='N']
        PE_size.append((abs(int(read.split('\t')[8])),max(obs_pos)-min(obs_pos))) # first is bowtie length, second is the length after all the trimming

    elif read_type=='UP':
        name,align = NGS_tools.parse_SAM(read,reference)
        alignment = NGS_tools.remove_lowQ_indels(align,Qcut)
        hiQ_indel = len([p for p in alignment if '-' in p])>0 
        if hiQ_indel: continue # if there is a high quality indel, skip and go to the next read


    # by this point there are no indels (first column of alignment==reference)
    alignment = [alignment[p] for p in coding] # trim to only the coding sequence (in frame)
    seqread = ''.join([p[1] if p[2]>=Qcut else 'N' for p in alignment])
    codons = sequence_tools.split_codons(seqread)

    # find the first and last codons that match the WT. Only use the codons between these. (this helps to minimize reads into the adaptor sequences)
    matches = [i for i in seqitr if codons[i]==WTcodons[i]]
    if len(matches)==0: continue
    codons = [i for i in enumerate(codons) if i[1] in code and (i[0]>=matches[0] and i[0]<=matches[-1])] # use codons that code for a specific amino acid and are within the valid range
    
    for i in codons:
        CDN_count[i[0]][i[1]]+=1
        AA_count[i[0]][code[i[1]]]+=1
       # for j in codons:
        #    if i[0]<j[0]:
         #       pair = code[i[1]]+code[j[1]]
          #      if pair_count[(i[0],j[0])].has_key(pair):
           #         pair_count[(i[0],j[0])][pair]+=1
            #    else:
             #       pair_count[(i[0],j[0])][pair]=1                    


outfile = filename.split('/')[-1][:-4]+'.pkl'
pickle.dump((CDN_count,AA_count,pair_count,PE_size),open(outfile,'wb'))
