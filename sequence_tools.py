from subprocess import Popen,STDOUT,PIPE
from random import choice
from os import remove
from time import time

## set the paths for various sequence utils
netblast_path = '/home/romero/code/netblast/netblast-2.2.21/bin/blastcl3'
clustalw_path = '/home/romero/code/clustalW/clustalw2'
muscle_path =   '/home/promero/code/muscle3.8.31/muscle3.8.31_i86linux64'
probcons_path = '/home/romero/code/probcons/probcons/probcons'
primer_path = '/home/romero/code/primer3-2.2.3/src/oligotm'
blastp_path = '/home/promero/code/ncbi-blast-2.2.28+/bin/blastp'

code = {'AAA': 'K','AAC': 'N','AAG': 'K','AAT': 'N','ACA': 'T','ACC': 'T','ACG': 'T','ACT': 'T','AGA': 'R','AGC': 'S','AGG': 'R','AGT': 'S','ATA': 'I','ATC': 'I','ATG': 'M','ATT': 'I',
        'CAA': 'Q','CAC': 'H','CAG': 'Q','CAT': 'H','CCA': 'P','CCC': 'P','CCG': 'P','CCT': 'P','CGA': 'R','CGC': 'R','CGG': 'R','CGT': 'R','CTA': 'L','CTC': 'L','CTG': 'L','CTT': 'L',
        'GAA': 'E','GAC': 'D','GAG': 'E','GAT': 'D','GCA': 'A','GCC': 'A','GCG': 'A','GCT': 'A','GGA': 'G','GGC': 'G','GGG': 'G','GGT': 'G','GTA': 'V','GTC': 'V','GTG': 'V','GTT': 'V',
        'TAA': '*','TAC': 'Y','TAG': '*','TAT': 'Y','TCA': 'S','TCC': 'S','TCG': 'S','TCT': 'S','TGA': '*','TGC': 'C','TGG': 'W','TGT': 'C','TTA': 'L','TTC': 'F','TTG': 'L','TTT': 'F',}


def rand_tag():
    '''creates a unique tag that can be used to identify the current files a script is working with.  necessary for running the same program multiple times in the same directory'''
    alpha = 'abcdefghijklmnopqrstuvwxyz0123456789'
    tag = ''.join([choice(alpha) for i in range(15)])
    return tag


####### BLAST tools ##########

def BLAST_search(query_seq): # the NCBI BLAST CHANGED!!
    fasta_str = '>query_sequence\n'+query_seq
    tag = rand_tag() # make a rand tag to make unique files
    fasta_filename = 'BLAST_'+tag+'.fasta'
    blast_filename = 'BLAST_'+tag+'.out'
    open(fasta_filename,'w').write(fasta_str)
    # run netBLAST
    BLAST_command = [blastp_path,
                     '-query '+fasta_filename, # input file
                     '-evalue 0.001', # set evalue.  10 is default (captures all)
                     '-db nr', # non-redundant database
                     '-max_target_seqs 100000', # don't want to be limited here
                     '-outfmt "10 mismatch sacc sseq"', # output as CSV and report: # of mismatches, the seq's GI, and the aligned portion's sequence
                     '-remote', # run on NCBI servers
                     '-out '+blast_filename] # set output file

    cmd = ' '.join(BLAST_command) # not sure why I have to do this?
    print(cmd)
    print ('submitting BLAST job')
    st = time()
    Popen(cmd,shell=True).wait() #runs netBLAST
    print ('BLAST search complete after %0.2f seconds' % (time()-st))
    remove(fasta_filename)
    return blast_filename


def GIlookup(GI):
    cmd = '/home/promero/code/edirect/edirect/efetch -db protein -id %s -format fasta' % str(GI) # this works for GI or Accession numbers
    output = Popen(cmd,shell=True,stdout=PIPE).communicate()[0]
    seq = ''.join(output.split('\n')[1:])
    return seq


def read_BLAST_output(filename):
    "requires a specific form of output file"
    data = open(filename,'r').read().split('\n\n\n\n\n')[1].split('Lambda')[0].strip().split('\n')
    data = [l for l in data if len(l)>0]
    seq_names = []
    for l in data:
        name = l.split()[0].strip()
        if name not in seq_names:
            seq_names.append(name)

    sequences = []
    for seq in seq_names:
        sequences.append(''.join([d.split()[2] for d in data if seq in d]))

    return seq_names,sequences



def read_BLAST_output_OLD(filename):
    "requires a specific form of output file"
    data = open(filename,'r').read().split('Query_1') # split on QUER, the first sequence of every section
    data = data[1:] # throw out header section
    data[-1] = data[-1].split('Lambda')[0] # throw out everything in the footer
    data = [d.strip() for d in data] # remove excess spaces and newlines
    data = [d.split('\n') for d in data] #split each section by \n
    data = [['FUCKDIS  '+l if i==0 else l for i,l in enumerate(d)] for d in data] # need this to get the spacing correct
    # take the middle section of each line
    d2 = []
    for ent in data:
        seqs = []
        for line in ent:
            blox = line.split()
            print(len(blox))
            if len(blox)==4:
                seqs.append(blox[2])
            elif len(blox)==2:
                seqs.append(blox[1])
            elif len(blox)==3:
                seqs.append(blox[2])
            else:
                seqs.append('')
                print(line)
        d2.append(seqs)

    data =  zip(*d2) # transpose the list of lists
    data = [''.join(d) for d in data] # join the segments
    #data = [d.replace('-','') for d in data] # remove all the gaps
    return data


# fucntions to run clustalW, muscle, probcons and read their output
def read_MSA_file(filename):
    "clustalW, muscle, and procons output fasta files"
    out_fasta = open(filename,'r').read()
    sequences = [''.join(seq.split('\n')[1:]) for seq in out_fasta.split('>')[1:]]
    alignment = zip(*sequences)
    return alignment


def clustalW_align(sequence_list,save=False,verbose=False):
    fasta_str = ''
    for i,seq in enumerate(sequence_list):
        fasta_str+='>seq'+str(i+1)+'\n'+seq+'\n\n'
    tag = rand_tag() # make a rand tag to make uniqe files
    fasta_filename = 'clustal_align_'+tag+'.fasta'
    clustal_filename = 'clustal_align_'+tag+'.out'
    open(fasta_filename,'w').write(fasta_str)
    clustal_command = [clustalw_path,
                       '-infile='+fasta_filename,
                       '-outorder=input', # keep sequences in order
                       '-outfile='+clustal_filename,
                       '-output=fasta']
    cmd = ' '.join(clustal_command) 
    if verbose:
        Popen(cmd,shell=True).wait() 
    else:
        Popen(cmd,shell=True,stderr=STDOUT, stdout=PIPE).wait()     
    remove(fasta_filename)
    remove('clustal_align_'+tag+'.dnd')
    alignment = read_MSA_file(clustal_filename)
    if save:
        return alignment,clustal_filename
    else:
        remove(clustal_filename)
        return alignment
    
    
def muscle_align(sequence_list,save=False,verbose=False):
    fasta_str = ''
    for i,seq in enumerate(sequence_list):
        fasta_str+='>seq'+str(i+1)+'\n'+seq+'\n\n'
    tag = rand_tag() # make a rand tag to make uniqe files
    fasta_filename = 'muscle_align_'+tag+'.fasta'
    temp_filename = 'muscle_align_'+tag+'.tmp'
    muscle_filename = 'muscle_align_'+tag+'.out'
    open(fasta_filename,'w').write(fasta_str)
    muscle_command = [muscle_path,
                      '-in '+fasta_filename, # input file
                      '-out '+temp_filename] # output file 
    cmd = ' '.join(muscle_command) # not sure why I have to do this?
    if verbose:
        Popen(cmd,shell=True).wait() 
    else:
        Popen(cmd,shell=True,stderr=STDOUT, stdout=PIPE).wait()     

    stable_command = ['python /home/promero/code/muscle3.8.31/stable.py',
                      fasta_filename, # input file
                      temp_filename, # alignment file
                      '>',
                      muscle_filename] # output file 

    Popen(' '.join(stable_command),shell=True,stderr=STDOUT, stdout=PIPE).wait()     

    remove(fasta_filename)
    remove(temp_filename)

    alignment = read_MSA_file(muscle_filename)    
    if save:
        return alignment,muscle_filename
    else:
        remove(muscle_filename)
        return alignment


def muscle_align_fast(sequence_list,save=False,verbose=False):
    fasta_str = ''
    for i,seq in enumerate(sequence_list):
        fasta_str+='>seq'+str(i+1)+'\n'+seq+'\n\n'
    tag = rand_tag() # make a rand tag to make uniqe files
    fasta_filename = 'muscle_align_'+tag+'.fasta'
    temp_filename = 'muscle_align_'+tag+'.tmp'
    muscle_filename = 'muscle_align_'+tag+'.out'
    open(fasta_filename,'w').write(fasta_str)
    muscle_command = [muscle_path,
                      '-in '+fasta_filename, # input file
                      '-out '+temp_filename, # output file
                      '-maxiters 1', '-diags -sv', '-distance1 kbit20_3'] # options for fast alignment
    cmd = ' '.join(muscle_command) # not sure why I have to do this?
    if verbose:
        Popen(cmd,shell=True).wait() 
    else:
        Popen(cmd,shell=True,stderr=STDOUT, stdout=PIPE).wait()     

    stable_command = ['python /home/promero/code/muscle3.8.31/stable.py',
                      fasta_filename, # input file
                      temp_filename, # alignment file
                      '>',
                      muscle_filename] # output file 

    Popen(' '.join(stable_command),shell=True,stderr=STDOUT, stdout=PIPE).wait()     

    remove(fasta_filename)
    remove(temp_filename)

    alignment = read_MSA_file(muscle_filename)    
    if save:
        return alignment,muscle_filename
    else:
        remove(muscle_filename)
        return alignment


def muscle_add_sequence(alignment,sequence,verbose=False):
    tag = rand_tag() # make a rand tag to make uniqe files
    new_seq = 'muscle_new_seq_'+tag+'.fasta'
    if isinstance(sequence,str): open(new_seq,'w').write('>new_seq\n'+sequence+'\n\n') ## just an ordinary sequence
    if isinstance(sequence,list) or isinstance(sequence,tuple): open(new_seq,'w').write(''.join(['>seq%i\n%s\n\n' % (i,''.join(seq)) for i,seq in enumerate(zip(*sequence))]))  ## in this case, adding an alignment
    existing_align = 'muscle_existing_align_'+tag+'.fasta'
    sequences = [''.join(s) for s in zip(*alignment)]
    open(existing_align,'w').write(''.join(['>seq'+str(i+1)+'\n'+sequences[i]+'\n\n' for i in range(len(sequences))]))
    muscle_filename = 'muscle_align_'+tag+'.out'
    muscle_command = [muscle_path,
                      '-profile ',
                      '-in1 '+existing_align,
                      '-in2 '+new_seq,
                      '-out '+muscle_filename] # output file 
    cmd = ' '.join(muscle_command) # not sure why I have to do this?
    if verbose:
        Popen(cmd,shell=True).wait() 
    else:
        Popen(cmd,shell=True,stderr=STDOUT, stdout=PIPE).wait()     
    alignment = read_MSA_file(muscle_filename)
    remove(muscle_filename)
    remove(new_seq)
    remove(existing_align)
    return alignment


def probcons_align(sequence_list,save=False):
    fasta_str = ''
    for i,seq in enumerate(sequence_list):
        fasta_str+='>seq'+str(i+1)+'\n'+seq+'\n\n'
    tag = rand_tag() # make a rand tag to make uniqe files
    fasta_filename = 'probcons_align_'+tag+'.fasta'
    probcons_filename = 'probcons_align_'+tag+'.out'
    open(fasta_filename,'w').write(fasta_str)
    probcons_command = [probcons_path,
                      fasta_filename, # input file
                      '> ' + probcons_filename] # output file 
    cmd = ' '.join(probcons_command) # not sure why I have to do this?
    Popen(cmd,shell=True).wait() #runs netBLAST
    remove(fasta_filename)
    alignment = read_MSA_file(probcons_filename)    
    if save:
        return alignment,probcons_filename
    else:
        remove(probcons_filename)
        return alignment


## fucntions to read, write, and print alignments
def read_alignment(filename):
    file = open(filename).read()
    data = [line for line in file.split('\n') if len(line) > 0 and line[0]!='#']
    if '>seq_names' in file:
        seq_names = data[data.index('>seq_names')+1:data.index('>alignment')]
    else:
        seq_names = []
    ali_data = data[data.index('>alignment')+1:]
    alignment =  [pos.split()[1:] for pos in ali_data]
    return alignment,seq_names


def write_alignment(alignment,seq_names=[],filename='new_alignment.aln'):
    num_digits = len(str(len(alignment)))
    alignment_file = open(filename,'w')
    if len(seq_names)>0:
        alignment_file.write('>seq_names\n')
        for name in seq_names:
            alignment_file.write(name+'\n')
    alignment_file.write('>alignment\n')
    for i,pos in enumerate(alignment):
        line = str(i).ljust(num_digits)+'  '+'  '.join(pos)
        alignment_file.write(line+'\n')
    alignment_file.close()


def print_alignment(alignment,seq_names=[]):
    if seq_names == []:
        seq_names=len(alignment[0])*['']
    name_length =  max([len(name) for name in seq_names])
    screen_width = 200
    num_lines = len(alignment)/screen_width    
    for i in range(num_lines+1):
        align_seg = alignment[(screen_width*i):(screen_width*(i+1))]
        conservation = []
        for pos in align_seg:
            if all([s==pos[0] for s in pos]):
                conservation.append('*')
            else:
                conservation.append(' ')
        seqs = zip(*align_seg)
        for i,seq in enumerate(seq_names):
            print (seq.ljust(name_length)+':'+''.join(seqs[i]))
        print (''.ljust(name_length)+':'+''.join(conservation)+'\n\n')


def quick_align(seq_list):
    align = muscle_align(seq_list)
    print_alignment(align)


def number_alignment(alignment):
    '''This numbers an alignment so the alignment position can be mapped back to the original sequence (and its features) that were aligned'''
    numbered_sequences = []
    for sequence in zip(*alignment):
        num_seq = []
        k = 0
        for pos in sequence:
            if pos=='-':
                num_seq.append('-')
            else:
                num_seq.append(k)
                k+=1
        numbered_sequences.append(num_seq)
    numbered_alignment = zip(*numbered_sequences)
    return numbered_alignment


def read_fasta(filename):
    data =  '\n'.join([l.strip() for l in open(filename).read().strip().split('\n') if len(l)>0 and l.strip()[0]!='#'])
    data = [s.strip() for s in data.split('>') if len(s)>0]
    seq_names = [d.split('\n')[0] for d in data]
    sequences = [''.join(d.split('\n')[1:]) for d in data]
    return seq_names,sequences


def read_fasta_dict(filename):
    data =  '\n'.join([l.strip() for l in open(filename).read().strip().split('\n') if len(l)>0 and l.strip()[0]!='#'])
    data = [s.strip() for s in data.split('>') if len(s)>0]
    sequences = dict((d.split('\n')[0],''.join(d.split('\n')[1:])) for d in data)
    return sequences



def write_fasta(seq_names,sequences,filename):
    open(filename,'w').write(''.join(['>%s\n%s\n\n' % (seq_names[i],sequences[i]) for i in range(len(sequences))]))



def read_clustal_output(filename): # for those people that run it online and want me to use it
    file = open(filename).read().replace('\r','')
    names = [l.split()[0] for l in file.split('\n\n')[0].split('\n')[:-1]]
    seqs = []
    for name in names:
        seqs.append(''.join([l.split()[1] for l in file.split('\n') if len(l)>0 and l.split()[0]==name]))
    alignment = zip(*seqs)
    return alignment,names


def translate(NT_seq):
    codons = [NT_seq[(3*i):(3*i)+3] for i in range(int(len(NT_seq)/3 ))]
#    AA_seq = ''.join(code[c] if code.has_key(c) else 'X' for c in codons)
    AA_seq = ''.join(code[c] if c in code else 'X' for c in codons)

    return AA_seq

def reverse_complement(NT_seq):
    num_seq = NT_seq.replace('A','1').replace('C','2').replace('G','3').replace('T','4')
    comp = list(num_seq.replace('1','T').replace('2','G').replace('3','C').replace('4','A'))
    comp.reverse()
    return ''.join(comp)


def est_Tm(NT_seq,mm=0):
    Tm = 81.5 + 41*float(NT_seq.count('G')+NT_seq.count('C'))/len(NT_seq) - (675.0/len(NT_seq)) - (100*float(mm)/len(NT_seq))
    return Tm



def oligo_Tm(primer):
    return float(Popen(primer_path+' '+primer,shell=True,stdout=PIPE).communicate()[0].strip())


def split_codons(NT_seq):
    codons = [NT_seq[(3*i):(3*i)+3] for i in range(len(NT_seq)//3 )]
    return codons


def assemble_sequences(seq1,seq2):
    """for assembling sequencing data"""
    seq1 = seq1+'XXX'
    ol_size = 5
    overlap = {}
    overlap[ol_size-1] = range(len(seq1))
    while len(overlap[ol_size-1])>1:
        overlap[ol_size] = []
        for i in overlap[ol_size-1]:
            if seq1[i:i+ol_size] in seq2:
                overlap[ol_size].append(i)
        ol_size += 1

    overlap_seq = seq1[overlap[ol_size-1][0]:overlap[ol_size-1][0]+ol_size-1] # this sequence is in both
    seq1 = seq1[:-3]
    assembly = []
    if seq1.find(overlap_seq) > seq2.find(overlap_seq):
        assembly.append(seq1[:seq1.find(overlap_seq)])    #put seq1 before
    else:
        assembly.append(seq2[:seq2.find(overlap_seq)])     #put seq2 before

    if len(seq1[seq1.find(overlap_seq):])>len(seq2[seq2.find(overlap_seq):]):
        assembly.append(seq1[seq1.find(overlap_seq):])
    else:
        assembly.append(seq2[seq2.find(overlap_seq):])

    sequence = ''.join(assembly)
    return sequence


def hamming_dist(s1,s2):
    """assumes s1 and s2 are the same length and aligned"""
    hd = len([i for i in range(len(s1)) if s1[i]!=s2[i]])
    return hd


def random_mutation(parent,num,alpha='ACEDGFIHKMLNQPSRTWVY'):
    mutant = ''.join([p for p in parent])
    while hamming_dist(parent,mutant) < num:
        pos = choice(range(len(mutant)))
        mutant = mutant[:pos]+choice(alpha)+mutant[pos+1:]
    return mutant



def print_translation(NT_seq):
    AA_seq = translate(NT_seq)
    print_alignment(zip(*(''.join([' '+p+' ' for p in AA_seq]),NT_seq)))


def pairwise_identity(alignment):
    sequences = zip(*alignment)
    identity = []
    for s1 in sequences:
        id = []
        for s2 in sequences:
            id.append(1 - float(hamming_dist(s1,s2))/len(s1))
        identity.append(id)
    return identity
            

def make_degenerate_code(cutoff=1e-6):
    IUBbases = {'A': {'A':1.0},'C': {'C':1.0},'G': {'G':1.0},'T': {'T':1.0},'R': {'A':0.5,'G':0.5},'Y': {'C':0.5,'T':0.5},'M': {'A':0.5,'C':0.5},'K': {'G':0.5,'T':0.5},'S': {'C':0.5,'G':0.5},'W': {'A':0.5,'T':0.5},'H': {'A':1.0/3,'C':1.0/3,'T':1.0/3},'B': {'C':1.0/3,'G':1.0/3,'T':1.0/3},'V': {'A':1.0/3,'C':1.0/3,'G':1.0/3},'D': {'A':1.0/3,'G':1.0/3,'T':1.0/3},'N': {'A':0.25,'C':0.25,'G':0.25,'T':0.25}}

    deg_codons = [b1+b2+b3 for b1 in IUBbases for b2 in IUBbases for b3 in IUBbases]

    deg_code = {}
    for codon in deg_codons:
        freq = {}
        for nt1 in IUBbases[codon[0]]:
            for nt2 in IUBbases[codon[1]]:
                for nt3 in IUBbases[codon[2]]:
                    AA = code[nt1+nt2+nt3]
                    prob = IUBbases[codon[0]][nt1]*IUBbases[codon[1]][nt2]*IUBbases[codon[2]][nt3]
                    if freq.has_key(AA):
                        freq[AA] += prob
                    else:
                        freq[AA] = prob
        deg_code[codon] = freq

    ## remove codons that have greater than 90% the same AA
    filtered = {}
    for codon in deg_code:
        if max(deg_code[codon].values())<0.90:
            filtered[codon] = deg_code[codon]

    deg_code = filtered 

    ## remove codons that have identical AA distributions
    filtered = {}
    observed = set()
    AA = sorted(set(code.values()))

    #choose the first codon just to get observed filled
    codon =deg_code.keys()[0]
    dist = tuple([deg_code[codon][a] if a in deg_code[codon] else 0 for a in AA])
    filtered[codon] = deg_code[codon]
    observed.add(dist)

    for codon in deg_code:
        dist = tuple([deg_code[codon][a] if a in deg_code[codon] else 0 for a in AA])
        if min([sum([(d[i]-dist[i])**2 for i in range(len(dist))])**0.5 for d in observed])>cutoff:
            filtered[codon] = deg_code[codon]
            observed.add(dist)

    return filtered



def designQCprimers(template,target):

    print( """from manual: 
          Primers should be between 25 and 45 bases in length, with a melting temperature (Tm) of >78 C.  Primers longer than 45 bases may be used, but using longer primers increases the likelihood of secondary structure
          The desired mutation (deletion or insertion) should be in the middle of the primer with 10-15 bases of correct sequence on both sides.
          The primers optimally should have a minimum GC content of 40% and should terminate in one or more C or G bases.""")


    align = muscle_align([template,target])
    matches = [i for i in range(len(align)) if align[i][0]==align[i][1]]
    mismatches = [i for i in range(len(align)) if i>min(matches) and i<max(matches) and align[i][0]!=align[i][1]] # these are the mismatches between the two matching ends

    primerrange = range(min(mismatches)-30,max(mismatches)+30)

    recentered = ''.join([align[i][1] if align[i][1]!='-' else align[i][0] for i in primerrange])
    align = muscle_align([template,recentered])
    matches = [i for i in range(len(align)) if align[i][0]==align[i][1]]
    mismatches = [i for i in range(len(align)) if i>min(matches) and i<max(matches) and align[i][0]!=align[i][1]] # these are the mismatches between the two matching ends

    primers = []
    for i in matches:
        for j in matches:
            if i<(j-20):
                primer = ''.join([a[1] for a in align[i:j]])

                # calc NT before and after
                before = min(mismatches)-i # number of bases before the first mismatch
                after = j-max(mismatches)-1 # number of bases after last mismatch
                unbalance = abs(before-after)

                GCend = primer[0] in 'GC' and primer[-1] in 'GC' # ends in GC
                goodlen = len(primer)>25 and len(primer)<45

                # calc Tm
                pGC = 100*float(primer.count('G')+primer.count('C'))/len(primer)
                pMM = 100*float(len(mismatches))/len(primer)
                Tm = 81.5 + 0.41*pGC - 675.0/len(primer) - pMM

                if pGC>40 and Tm>78 and Tm<90 and GCend and goodlen and unbalance<5:
                    primers.append((unbalance,Tm,pGC,len(primer),i,j,before,after))

    primers = sorted(primers)
    for k,pr in enumerate(primers):
        unbalance,Tm,pGC,plen,i,j,before,after = pr
        print ('PRIMER #',k+1)
        print ('length: %i, before/after: %i/%i, Tm: %0.1f, percent GC: %0.1f'%(plen,before,after, Tm, pGC))
        print_alignment(align[i:j])


def align2phylip(alignment,names=''):
    sequences = [''.join(s) for s in zip(*alignment)]
    if names=='': names = ['seq%i'%(i+1) for i in range(len(sequences))]
    if len(names)!=len(sequences): print ('problem'),1/0

    numseq = len(sequences)
    numpos = len(alignment)
    maxnamelen = max([len(n) for n in names])

    phylip = ' %i %i \n' % (numseq,numpos)
    phylip += '\n'.join([names[i].ljust(maxnamelen+10)+sequences[i] for i in range(numseq)])

    return phylip



#From: http://www.genscript.com/cgi-bin/tools/codon_freq_table
#Fields: Triplet | Amino acid | Fraction | Frequency/Thousand | (Number)
cdn_usage = """TTT F 0.58 22.1( 80995)  TCT S 0.17 10.4( 38027)  TAT Y 0.59 17.5( 63937)  TGT C 0.46  5.2( 19138)
TTC F 0.42 16.0( 58774)  TCC S 0.15  9.1( 33430)  TAC Y 0.41 12.2( 44631)  TGC C 0.54  6.1( 22188)
TTA L 0.14 14.3( 52382)  TCA S 0.14  8.9( 32715)  TAA * 0.61  2.0(  7356)  TGA * 0.30  1.0(  3623)
TTG L 0.13 13.0( 47500)  TCG S 0.14  8.5( 31146)  TAG * 0.09  0.3(   989)  TGG W 1.00 13.9( 50991)

CTT L 0.12 11.9( 43449)  CCT P 0.18  7.5( 27340)  CAT H 0.57 12.5( 45879)  CGT R 0.36 20.0( 73197)
CTC L 0.10 10.2( 37347)  CCC P 0.13  5.4( 19666)  CAC H 0.43  9.3( 34078)  CGC R 0.36 19.7( 72212)
CTA L 0.04  4.2( 15409)  CCA P 0.20  8.6( 31534)  CAA Q 0.34 14.6( 53394)  CGA R 0.07  3.8( 13844)
CTG L 0.47 48.4(177210)  CCG P 0.49 20.9( 76644)  CAG Q 0.66 28.4(104171)  CGG R 0.11  5.9( 21552)

ATT I 0.49 29.8(109072)  ACT T 0.19 10.3( 37842)  AAT N 0.49 20.6( 75436)  AGT S 0.16  9.9( 36097)
ATC I 0.39 23.7( 86796)  ACC T 0.40 22.0( 80547)  AAC N 0.51 21.4( 78443)  AGC S 0.25 15.2( 55551)
ATA I 0.11  6.8( 24984)  ACA T 0.17  9.3( 33910)  AAA K 0.74 35.3(129137)  AGA R 0.07  3.6( 13152)
ATG M 1.00 26.4( 96695)  ACG T 0.25 13.7( 50269)  AAG K 0.26 12.4( 45459)  AGG R 0.04  2.1(  7607)

GTT V 0.28 19.8( 72584)  GCT A 0.18 17.1( 62479)  GAT D 0.63 32.7(119939)  GGT G 0.35 25.5( 93325)
GTC V 0.20 14.3( 52439)  GCC A 0.26 24.2( 88721)  GAC D 0.37 19.2( 70394)  GGC G 0.37 27.1( 99390)
GTA V 0.17 11.6( 42420)  GCA A 0.23 21.2( 77547)  GAA E 0.68 39.1(143353)  GGA G 0.13  9.5( 34799)
GTG V 0.35 24.4( 89265)  GCG A 0.33 30.1(110308)  GAG E 0.32 18.7( 68609)  GGG G 0.15 11.3( 41277)"""
cdn_usage = dict((c,float(cdn_usage[cdn_usage.find(c)+6:cdn_usage.find(c)+10])) for c in code)
#####
import re
from Bio import SeqIO
def find_max_orf(fasta_file):

    seq = SeqIO.read(fasta_file, 'fasta').seq #load sequnece from fasta
    seq = str(seq).upper() #convert to string, capitalize entire sequence
    
    orf = max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',seq), key = len)
    start = re.search(orf,seq).start()
    stop = re.search(orf,seq).end()
    #   find longest ORF
    print(len(orf))
    print(start, stop)
    return start, stop, orf