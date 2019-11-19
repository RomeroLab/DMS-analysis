



# load read counts
read_counts = open('Bgl3_DMS_first1000reads_read_counts.txt').read().strip().split('\n')

# flatten mutation file into a single list of mutations 
mutations = open('Bgl3_DMS_first1000reads_mutations.txt').read().strip().replace('\n',',').split(',')


AAs = 'ACDEFGHIKLMNPQRSTVWY*'


counts = []
for i in range(len(read_counts)):
    ind,AA,rc,mc = read_counts[i].replace(' ','').split(',')
    WT_count = int(rc) - int(mc)
    cnt = []
    for a in AAs:
        if a==AA:
            cnt.append(WT_count)
        else:
            mut = AA+str(i)+a
            cnt.append(mutations.count(mut))
    counts.append(cnt)
