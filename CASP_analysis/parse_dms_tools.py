#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 18:02:44 2018

@author: hridindu
"""

import pickle
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#import matplotlib.pyplot as plt
import numpy as np
import sequence_tools as st
import pandas as pd
#import scipy as sp
#import math
#import heapq
#from scipy import stats
np.seterr(divide = 'ignore') 
#%%

def parse_DMS(infile, ref_seq_fasta, coding_start, coding_term): #Input .pkl file as infile, reference fasta sequence as refseq
    
#initialize the reference sequence. 
#This is probably the fasta file of the insert that you sent for sequnceing, 
#and wil likely have multiple reading frames. Make sure to specify the correcty reading frame (0-indexed).  
    
    ref_seq_handle = open(ref_seq_fasta, 'rU')
    ref_rec = list(SeqIO.parse(ref_seq_handle, "fasta"))
    ref_seq = str(ref_rec[0].seq)
    coding = ref_seq[coding_start:coding_term]
    aa_seq = st.translate(coding)
    codon_count = int(len(coding)/3)
    ref_codons =[]
    for n in range(codon_count):
        ref_codons.append(coding[n*3:(n*3)+3])
    
    
    CDN_count,AA_count,pair_count,PE_size = pickle.load (open(infile, 'rb'))
    return CDN_count,AA_count,pair_count,PE_size, ref_codons, aa_seq
#%%    
def read_cov_DMS(infile, ref_seq_fasta, coding_start, coding_term):
#    with open(infile, 'rb') as infile:
    CDN_count,AA_count,pair_count,PE_size,ref_codons, aa_seq = parse_DMS(infile, ref_seq_fasta, coding_start, coding_term)
    
    read_coverage = []
    NT_mut_freq = []    
#    trimmed = [sizes[1] for sizes in PE_size]
#    untrimmed = [sizes[0] for sizes in PE_size]
    
    for pos in CDN_count:
        total = 0
        mut = 0
        for cdn in CDN_count[pos]:
            if 'N' not in cdn:
                total += CDN_count[pos][cdn]
                if cdn != ref_codons[pos]:
                    mut += CDN_count[pos][cdn] 
        NT_mut_freq.append(np.log2(mut/total))
        read_coverage.append(total /1000)
    return read_coverage, NT_mut_freq
#    infile.close()
 #%%   
def aa_data_DMS(infiles, ref_seq_fasta, coding_start, coding_term, offset):
    
    aa_count_data = []
    aa_freq_data = []
    post_mut_freq_data = []
    
    for infile in infiles:
        
        CDN_count,AA_count,pair_count,PE_size,ref_codons, aa_seq = parse_DMS(infile, ref_seq_fasta, coding_start, coding_term)
        reads,_ = read_cov_DMS(infile, ref_seq_fasta, coding_start, coding_term)
        def calc_freq(AA_count,thresh=10):
            AAs = sorted(AA_count[0].keys())
            count = np.array([[AA_count[p][a] for a in AAs] for p in AA_count])
            count[count<thresh] = 0
            freq = count/np.tile(count.sum(1),(21,1)).T
            return freq
       
        flat_freq = [float('NaN') if (sum(AA_count[i].values()) == 0 or AA_count[i].values() == 0)
                    else (np.log2((1 - max(AA_count[i].values())/sum(AA_count[i].values())))) for i in AA_count]                  
        
        pos_mut_freq = pd.DataFrame(data = flat_freq, 
                                    index = range(offset, len(aa_seq) + offset))
        aa_count = pd.DataFrame(data = AA_count, 
                                columns = range(offset, len(aa_seq) + offset)).T
        aa_freq = pd.DataFrame(calc_freq(AA_count), 
                               columns = list(aa_count.columns.values),
                               index = range(offset, len(aa_seq) + offset))
#        
#        pos_mut_freq = pd.DataFrame(data = flat_freq, 
#                                    index = range( len(aa_seq)))
#        aa_count = pd.DataFrame(data = AA_count, 
#                                columns = range( len(aa_seq))).T
#        aa_freq = pd.DataFrame(calc_freq(AA_count), 
#                               columns = list(aa_count.columns.values),
#                               index = range(len(aa_seq)))
        
        aa_count['WT'] = list(aa_seq)
        aa_freq['WT'] = list(aa_seq)
        pos_mut_freq['WT'] = list(aa_seq)
        
        if 'PS' in infile:
            aa_count['Screened'] = False
            aa_freq['Screened'] = False
            pos_mut_freq['Screened'] = False
        else: 
            aa_count['Screened'] = True
            aa_freq['Screened'] = True
            pos_mut_freq['Screened'] = True
            
        aa_count['Date'] = infile.split('.')[0].split('P')[0]
        aa_freq['Date'] = infile.split('.')[0].split('P')[0]
        pos_mut_freq['Date'] = infile.split('.')[0].split('P')[0]
        
        aa_count['Gene'] = open(ref_seq_fasta, 'r').readline().strip()
        aa_freq['Gene'] = open(ref_seq_fasta, 'r').readline().strip()
        pos_mut_freq['Gene'] = open(ref_seq_fasta, 'r').readline().strip()
        
        aa_count['Position'] = list(range(offset, len(aa_seq) + offset))
        aa_freq['Position'] = list(range(offset, len(aa_seq) + offset))
        pos_mut_freq['Position'] = list(range(offset, len(aa_seq) + offset))
        
        aa_freq['read_coverage'] = reads
                
        aa_count_data.append(aa_count)
        aa_freq_data.append(aa_freq)
        post_mut_freq_data.append(pos_mut_freq)
            
    aa_count = pd.concat(aa_count_data)
    aa_freq = pd.concat(aa_freq_data)
    pos_mut_freq = pd.concat(post_mut_freq_data)
    
    res = []

    for infile in infiles:
        if 'PS' in infile:
            pass
        else:
            CDN_count,AA_count,pair_count,PE_size,ref_codons, aa_seq = parse_DMS(infile, ref_seq_fasta, coding_start, coding_term)
    
            Date = infile[:6]
            pri = aa_freq.loc[(aa_freq['Date'] == Date) & (aa_freq['Screened'] == False)].iloc[:,:21]
            post = aa_freq.loc[(aa_freq['Date'] == Date) & (aa_freq['Screened'] == True)].iloc[:,:21]
            
            pri = pri.select_dtypes(include = ['float64']).replace([np.inf, -np.inf, 0], np.nan)
            post = post.select_dtypes(include = ['float64']).replace([np.inf, -np.inf, 0], np.nan)
        
            re = pd.DataFrame(data = post * np.log2(post / pri))
            re_flat = pd.DataFrame({'RE' : list(re.sum(axis = 1))}, index = range(offset, len(aa_seq) + offset))
            
            re_flat = re.sum(axis = 1, skipna = True)
            
#            euclid = pd.DataFrame(data = (np.log2((post**2/pri**2)**(1/2))).sum(axis = 1))
            euclid = pd.DataFrame(data = (((post**2-pri**2)**(1/2))).sum(axis = 1, skipna = True))
            
            re['WT'] = list(aa_seq)
            re['Date'] = Date
            re['Gene'] = open(ref_seq_fasta, 'r').readline().strip()
#            re['RE_flat'] = (re_flat - re_flat.min()) / (re_flat.max() - re_flat.min())
            re['RE_flat'] = re_flat
            re['euclid'] = euclid
            re['pct_RE'] = re['RE_flat'].rank(pct = True)
            re['Position'] = re.index
#            re = re.replace([np.inf, -np.inf], np.nan).fillna(0)
            
            res.append(re)

    re = pd.concat(res).reset_index(drop = True)
    
    return aa_count, aa_freq, pos_mut_freq, re

#%%
def calc_DMS_RE(post_file, pre_file, ref_seq_fasta, coding_start, coding_term):
    
    aa_count_post, aa_freq_post, pos_mut_freq_post = aa_data_DMS(post_file, ref_seq_fasta, coding_start, coding_term)
    aa_count_pre, aa_freq_pre, pos_mut_freq_pre = aa_data_DMS(pre_file, ref_seq_fasta, coding_start, coding_term)
    
    re = aa_freq_post * np.log2(aa_freq_post/aa_freq_pre)
    nans = np.isnan(re)
    infs = np.isinf(re)
    re[infs] = 0
    re[nans] = 0
    flat_re = re.sum(axis = 1)
    return re, flat_re


#%%
def read_seq(ref_seq_fasta, coding_start, coding_term, offset):
    ref_seq_handle = open(ref_seq_fasta, 'rU')
    ref_rec = list(SeqIO.parse(ref_seq_handle, "fasta"))
    ref_seq = str(ref_rec[0].seq)
    coding = ref_seq[coding_start:coding_term]
    aa_seq = st.translate(coding)
    
    return aa_seq
#%%


from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62
gap_open = -10
gap_extend = -0.5

def pw_align(seq1, seq2): #input two amino acid sequences. You can use read_seq to generate them from reference fasta files. 
    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    top_aln = alns[0]

    aln1, aln2, score, begin, end = top_aln
    print(aln1+'\n'+aln2)

#%%

def write_mut_list(pickle_infile, text_outfile):
    
    with open(pickle_infile, 'rb') as mut_file:
        mutations = pickle.load(mut_file)

        outfile = open(str(text_outfile) + '.txt', 'w+')
        outfile.write('\n'.join(map(str, mutations)))
        mut_file.close()
    return text_outfile







