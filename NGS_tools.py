
# this is the Illumina Qscore mapping                                                                                                                                                       
Qscore = dict((chr(i),i-33) for i in range(33,90))


def parse_SAM(read,reference): # inputs a line from a SAM file and the reference sequence
    # organize data
    data = read.split('\t')
    name = data[0]
    start = int(data[3])-1 # get into python numbering
    cigar = data[5]
    seq = data[9]
    Q = [Qscore[q] for q in data[10]]

    #cigar2clist: split the cigar string into sections
    clist = []
    pos = ''
    for p in cigar:
        if p.isalpha():
            clist.append((int(pos),p))
            pos = ''
        else:
            pos+=p

    # clip the seq and Q according to the clist
    if clist[0][1]=='S':
        seq = seq[clist[0][0]:]
        Q = Q[clist[0][0]:]
        clist.pop(0)
    if clist[-1][1]=='S':
        seq = seq[:-clist[-1][0]]
        Q = Q[:-clist[-1][0]]
        clist.pop(-1)

    # cycle through clist and generate alignment
    refstart = start
    readstart = 0
    alignment = [(reference[i],'N',0) for i in range(start)] # first part of alignment is just reference and Ns
    for section in clist:
        if section[1]=='M': # match: add to both
            ref = reference[refstart:refstart+section[0]]
            read = seq[readstart:readstart+section[0]]
            readQ = Q[readstart:readstart+section[0]]

            alignment.extend(zip(ref,read,readQ))

            refstart += section[0]
            readstart += section[0]

        if section[1]=='D': # deletion: only add to ref seq
            ref = reference[refstart:refstart+section[0]]
            read = '-'*section[0]
            readQ = '-'*section[0]

            alignment.extend(zip(ref,read,readQ))

            refstart += section[0]

        if section[1]=='I': # insertion: only add to read seq
            ref = '-'*section[0]
            read = seq[readstart:readstart+section[0]]
            readQ = Q[readstart:readstart+section[0]]

            alignment.extend(zip(ref,read,readQ))

            readstart += section[0]

    alignment.extend([(reference[i],'N',0) for i in range(refstart,len(reference))]) # add on the remainder of the referenbce sequence

    # need to assign Qscores for the deletions: use the surrounding base Qscores and mean or min
    
    deletions = [i for i in range(len(alignment)) if alignment[i][2]=='-']
    for pos in deletions:
                
        Qbefore = alignment[pos-1][2]
#        Qafter = alignment[pos+1][2]
#        delQ = (Qbefore+Qafter)/2 # take the average of the Q socres before and after (could also take the min)
#        delQ = min([Qbefore,Qafter])
        delQ = 40
        alignment[pos] = alignment[pos][:2]+(delQ,)

    return name, alignment


def remove_lowQ_indels(alignment,Qcut):
    # remove low quality insertions
    alignment = [p for p in alignment if p[0]!='-' or (p[0]=='-' and p[2]>=Qcut)]#

    # set low quality deletions to WT
    alignment = [(p[0],p[0],p[2]) if (p[1]=='-' and p[2]<Qcut) else p for p in alignment]

    return alignment

