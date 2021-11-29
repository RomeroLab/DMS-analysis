#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 09:53:41 2019

@author: hridindu
"""

import matplotlib.pyplot as plt
import numpy as np
import sequence_tools as st
import pandas as pd
from scipy import stats

import parse_dms_tools as dms
import parse_msa_tools as msa
#import parse_dssp as dssp
import seaborn as sns; sns.set()
sns.set_style('white')
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO

my_path = '/Users/hridindu/Documents/caspase_project/121018_illumina/analysis_pass_4'

#reference sequences with subsequent coding start, coding terminus, n-terminus offset. 
c3 = 'caspase3_ref.fasta'
(c3o, c3t, c3_off) = (149,926,28)

c7 = 'caspase7_ref.fasta'
(c7o, c7t, c7_off) = (148, 994, 23)

c3_datasets =['092118PS.pkl', 
              '100518PS.pkl', 
              '100618PS.pkl',
              '092118.pkl', 
              '100518.pkl', 
              '100618.pkl', ]

c7_datasets =['101518PS.pkl', 
              '101618PS.pkl', 
              '102618PS.pkl',
              '101518.pkl', 
              '101618.pkl', 
              '102618.pkl' ]
#This list is easiest to iterate over for looking at relative entropies, just because you don't have to deal with
#redundant dataset names
c3_sorts = ['092118.pkl', 
            '100518.pkl', 
            '100618.pkl', ]

c7_sorts = ['101518.pkl', 
            '101618.pkl', 
            '102618.pkl' ]

c3_aa_count, c3_aa_freq, c3_pos_mut_freq, c3_re = dms.aa_data_DMS(c3_datasets,c3,c3o,c3t,c3_off)
c7_aa_count, c7_aa_freq, c7_pos_mut_freq, c7_re = dms.aa_data_DMS(c7_datasets,c7,c7o,c7t,c7_off)

c3_msa_aa, _ = msa.count_MSA('C3_uniref90_filtered.aln', c3_off)
c7_msa_aa, _ = msa.count_MSA('C7_uniref90_filtered.aln', c7_off)
c3_msa_re = msa.calc_MSA_RE('C3_uniref90_filtered.aln', c3_off)
c7_msa_re = msa.calc_MSA_RE('C7_uniref90_filtered.aln', c7_off)

"""Load the Caspase-3 Caspase-7 Sequence alignment"""
aln = list(SeqIO.parse('c3_c7_aln.fasta', 'fasta'))

c3_aln = {
        'WT' : [aa for pos, aa in enumerate(aln[1].seq + '*', 1) if aa != '-'],
        'aln_pos' : [pos for pos, aa in enumerate(aln[1].seq + '*', 1) if aa != '-'],
        'gene': [aln[1].id for pos, aa in enumerate(aln[1].seq + '*', 1) if aa != '-']
            }
c7_aln = {
        'WT' : [aa for pos, aa in enumerate(aln[0].seq + '*' , 1) if aa != '-'],
        'aln_pos' : [pos for pos, aa in enumerate(aln[0].seq + '*' , 1) if aa != '-'],
        'gene': [aln[0].id for pos, aa in enumerate(aln[0].seq + '*', 1) if aa != '-']
            }

c3_aln = pd.DataFrame(data = c3_aln, dtype = 'float')
c7_aln = pd.DataFrame(data = c7_aln, dtype = 'float')

c3_aln_long = pd.concat([c3_aln,c3_aln,c3_aln], axis = 0)
c7_aln_long = pd.concat([c7_aln,c7_aln,c7_aln], axis =0)

c3_re['Aligned Position'] = list(c3_aln_long['aln_pos'])
c7_re['Aligned Position'] = list(c7_aln_long['aln_pos'])

c3_mean = pd.DataFrame(data = {'Position' : c3_re['Position'].unique(),
                               'Aligned Position' : c3_aln['aln_pos'],
                               'WT Sequence' : c3_aln['WT'],
                               'mean_RE': c3_re.loc[c3_re['Date'] != '102116'].groupby('Position')['RE_flat'].mean().reset_index(drop = True),
                               'std_err_re' : c3_re.loc[c3_re['Date'] != '102116'].groupby('Position')['RE_flat'].sem().reset_index(drop = True),
                              })
                            

                                 
                        
c7_mean = pd.DataFrame(data = {'Position' : c7_re['Position'].unique(),
                               'Aligned Position' : c7_aln['aln_pos'],
                               'WT Sequence' : c7_aln['WT'],
                               'mean_RE': c7_re.groupby('Position')['RE_flat'].mean().reset_index(drop = True),
                               'std_err_re' : c7_re.groupby('Position')['RE_flat'].sem().reset_index(drop = True)
                              })
                                 

c3_mean = c3_mean.fillna(0)
c7_mean = c7_mean.fillna(0)

c_mean = pd.merge(c3_mean, c7_mean, how = 'outer',  on = 'Aligned Position', suffixes= ('_c3', '_c7'))
c_mean['shared_res_color'] = ['r' if aa == bb else 'b' for (aa,bb) in (zip(c_mean['WT Sequence_c3'],c_mean['WT Sequence_c7']))]

mw_pvals = []
t_pvals = []
for i, pos in enumerate(c_mean['Aligned Position'],1):
    c3 = list(c3_re.loc[(c3_re['Aligned Position'] == pos), 'RE_flat'])
    c7 = list(c7_re.loc[(c7_re['Aligned Position'] == pos), 'RE_flat'])
    if (c3 == c7 or len(c3) != len(c7) or np.isnan(c3).any() or np.isnan(c7).any()):
        mw_pvals.append(float('NaN'))
        t_pvals.append(float('NaN'))
    else:
        uval, mw_pval = stats.mannwhitneyu(c3,c7)
        mw_pvals.append(mw_pval)
        tval, t_pval = stats.ttest_ind(c3,c7, equal_var=False)
        t_pvals.append(t_pval)
c_mean['mw_p_value'] = mw_pvals
c_mean['t_p_value'] = t_pvals

#%%
c_mean['RE_deviation'] = np.log10(c_mean['mean_RE_c3'] / c_mean['mean_RE_c7'])

#%%
def three_to_one(three_letter_code):
    mapping = {'Aba':'A','Ace':'X','Acr':'X','Ala':'A','Aly':'K','Arg':'R','Asn':'N','Asp':'D','Cas':'C',
           'Ccs':'C','Cme':'C','Csd':'C','Cso':'C','Csx':'C','Cys':'C','Dal':'A','Dbb':'T','Dbu':'T',
           'Dha':'S','Gln':'Q','Glu':'E','Gly':'G','Glz':'G','His':'H','Hse':'S','Ile':'I','Leu':'L',
           'Llp':'K','Lys':'K','Men':'N','Met':'M','Mly':'K','Mse':'M','Nh2':'X','Nle':'L','Ocs':'C',
           'Pca':'E','Phe':'F','Pro':'P','Ptr':'Y','Sep':'S','Ser':'S','Thr':'T','Tih':'A','Tpo':'T',
           'Trp':'W','Tyr':'Y','Unk':'X','Val':'V','Ycm':'C','Sec':'U','Pyl':'O'} # you can add more
    return mapping[three_letter_code[0].upper() + three_letter_code[1:].lower()]



#%%
with open('2ql5.pdb', 'r') as infile,  open('2ql5_RE_dev.pdb' ,'w') as outfile:
        
        for line in infile:
                
            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_mean.loc[c_mean['Position_c7'] == int(pdb_pos), 'WT Sequence_c7']:
                    print('fuck up, aa seq doesnt align')
                
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c7'] == int(pdb_pos), 'RE_deviation'])).rjust(6)
                if 'nan' in re:
                    re = str('%0.2f' % 100).rjust(6)
                    line.replace(b,re)
                else:
                        line = line.replace(b, re)
                        outfile.write(line)
            
            if line[:4] == ('ATOM') and (line[21] in 'CD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                if aa in c_mean.loc[c_mean['Position_c7'] == int(pdb_pos)-300, 'WT Sequence_c7']:
                    print('fuck up, aa seq doesnt align')
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c7'] == int(pdb_pos)-300, 'RE_deviation'])).rjust(6)
                if 'nan' in re:
                    re = str('%0.2f' % 100).rjust(6)
                    line.replace(b,re)
                else:                        
                    line = line.replace(b, re)
                    outfile.write(line)
           
            else: 
                outfile.write(line)
outfile.close()
#%%

with open('2h5j.pdb', 'r') as infile,  open('2h5j_RE_dev.pdb' ,'w') as outfile:
        
        
        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 'WT Sequence_c3']:
                    print('fuck up, aa seq doesnt align')
                
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 'RE_deviation'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
            
            if line[:4] == ('ATOM') and (line[21] in 'CD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                if aa in c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 'WT Sequence_c3']:
                    print('fuck up, aa seq doesnt align')
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 'RE_deviation'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
           
            else: 
                outfile.write(line)
outfile.close()

#%%
with open('2ql5.pdb', 'r') as infile,  open('2ql5_RE.pdb' ,'w') as outfile:

        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_mean.loc[c_mean['Position_c7'] == int(pdb_pos), 'WT Sequence_c7']:
                    print('fuck up, aa seq doesnt align')
                
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c7'] == int(pdb_pos), 'mean_RE_c7'])).rjust(6)
                if 'nan' in re:
                    re = str('%0.2f' % 100).rjust(6)
                    line.replace(b,re)
                else:                        
                    line = line.replace(b, re)
                    outfile.write(line)
            
            if line[:4] == ('ATOM') and (line[21] in 'CD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                if aa in c_mean.loc[c_mean['Position_c7'] == int(pdb_pos)-300, 'WT Sequence_c7']:
                    print('fuck up, aa seq doesnt align')
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c7'] == int(pdb_pos)-300, 'mean_RE_c7'])).rjust(6)
                if 'nan' in re:
                    re = str('%0.2f' % 100).rjust(6)
                    line.replace(b,re)
                else:                        
                    line = line.replace(b, re)
                    outfile.write(line)
           
            else: 
                outfile.write(line)
outfile.close()
#%%
with open('2h5j.pdb', 'r') as infile,  open('2h5j_RE.pdb' ,'w') as outfile:

        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'ABCD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 'WT Sequence_c3']:
                    print('fuck up, aa seq doesnt align')
                
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 'mean_RE_c3'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
           
            else: 
                outfile.write(line)
outfile.close()
#%%

with open('2ql5.pdb', 'r') as infile,  open('2ql5_pval.pdb' ,'w') as outfile:
        
        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_mean.loc[c_mean['Position_c7'] == int(pdb_pos), 'WT Sequence_c7']:
                    print('fuck up, aa seq doesnt align')
                
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c7'] == int(pdb_pos), 't_p_value'])).rjust(6)
                if 'nan' in re:
                    re = str('%0.2f' % 100).rjust(6)
                    line.replace(b,re)
                else:                        
                    line = line.replace(b, re)
                    outfile.write(line)
            
            if line[:4] == ('ATOM') and (line[21] in 'CD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                if aa in c_mean.loc[c_mean['Position_c7'] == int(pdb_pos)-300, 'WT Sequence_c7']:
                    print('fuck up, aa seq doesnt align')
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c7'] == int(pdb_pos)-300, 't_p_value'])).rjust(6)
                if 'nan' in re:
                    re = str('%0.2f' % 100).rjust(6)
                    line.replace(b,re)
                else:                        
                    line = line.replace(b, re)
                    outfile.write(line)
            else: 
                outfile.write(line)
outfile.close()
#%%
with open('2h5j.pdb', 'r') as infile,  open('2h5j_pval.pdb' ,'w') as outfile:

        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'ABCD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 'WT Sequence_c3']:
                    print('fuck up, aa seq doesnt align')
                
                re = str('%0.2f' % np.abs(c_mean.loc[c_mean['Position_c3'] == int(pdb_pos), 't_p_value'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
           
            else: 
                outfile.write(line)
outfile.close()


#%%

















































