#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 10:54:55 2020

@author: hridindu
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import sequence_tools as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats

template = 'none'

import parse_dms_tools as dms
import parse_msa_tools as msa
#import parse_dssp as dssp
import seaborn as sns; sns.set()
sns.set_style('white')

import pickle
from sklearn import preprocessing
import scipy.stats
min_max_scaler = preprocessing.MinMaxScaler(feature_range=(1e-6, 100))

my_path = '/Users/hridindu/Documents/C3-C7_paper/PU'
#reference sequences with subsequent coding start, coding terminus, n-terminus offset. 
c3 = 'caspase3_ref.fasta'
(c3o, c3t, c3_off) = (149,926,28)

c7 = 'caspase7_ref.fasta'
(c7o, c7t, c7_off) = (148, 994, 23)

important_positions = [164, 122, 176]
div_mask = [173, 230, 163, 184, 220, 224]

def three_to_one(three_letter_code):
    mapping = {'Aba':'A','Ace':'X','Acr':'X','Ala':'A','Aly':'K','Arg':'R','Asn':'N','Asp':'D','Cas':'C',
           'Ccs':'C','Cme':'C','Csd':'C','Cso':'C','Csx':'C','Cys':'C','Dal':'A','Dbb':'T','Dbu':'T',
           'Dha':'S','Gln':'Q','Glu':'E','Gly':'G','Glz':'G','His':'H','Hse':'S','Ile':'I','Leu':'L',
           'Llp':'K','Lys':'K','Men':'N','Met':'M','Mly':'K','Mse':'M','Nh2':'X','Nle':'L','Ocs':'C',
           'Pca':'E','Phe':'F','Pro':'P','Ptr':'Y','Sep':'S','Ser':'S','Thr':'T','Tih':'A','Tpo':'T',
           'Trp':'W','Tyr':'Y','Unk':'X','Val':'V','Ycm':'C','Sec':'U','Pyl':'O'} # you can add more
    return mapping[three_letter_code[0].upper() + three_letter_code[1:].lower()]
#%%
c3 = pickle.load(open('/Users/hridindu/Documents/C3-C7_paper/PU/C3.pkl', 'rb'))
c7 = pickle.load(open('/Users/hridindu/Documents/C3-C7_paper/PU/C7.pkl', 'rb'))
c_mean = pickle.load(open('/Users/hridindu/Documents/C3-C7_paper/PU/agg_casp_PU.pkl', 'rb'))

mask = c_mean.loc[c_mean['Aligned Position'].isin(important_positions)]
#%%

sig = 1
c_mean = c_mean.loc[((c_mean['p.grp.adj_c3'] < sig) & (c_mean['p.grp.adj_c7'] < sig))]

####   Make df with mean(abs(coef))
# mean_abs_coef_c7 = c_mean[['coef_c7', 
#                            'Position_c7', 
#                            'WT_aa_c7', 
#                            'Aligned Position']].groupby(['Position_c7', 'WT_aa_c7', 'Aligned Position']).apply(lambda x: x.abs().mean()).reset_index()
# mean_abs_coef_c3 = c_mean[['coef_c3', 
#                            'Position_c3', 
#                            'WT_aa_c3', 
#                            'Aligned Position']].groupby(['Position_c3', 'WT_aa_c3', 'Aligned Position']).apply(lambda x: x.abs().mean()).reset_index()


mean_abs_coef_c7 = c_mean[['coef_c7', 
                           'Position_c7', 
                           'WT_aa_c7', 
                           'Aligned Position', 
                           'mean_RE_c7',
                           'p_c7']].groupby(['Position_c7', 'WT_aa_c7', 'Aligned Position','mean_RE_c7']).apply(lambda x: x.abs().mean()).reset_index()

mean_abs_coef_c3 = c_mean[['coef_c3', 
                            'Position_c3', 
                            'WT_aa_c3', 
                            'Aligned Position',
                            'mean_RE_c3',
                            'p_c3']].groupby(['Position_c3', 'WT_aa_c3', 'Aligned Position','mean_RE_c3']).apply(lambda x: x.abs().mean()).reset_index()

#%%
c_coef_div = mean_abs_coef_c3.merge(mean_abs_coef_c7, on = 'Aligned Position')
c_coef_div['div'] = c_coef_div['coef_c3'] - c_coef_div['coef_c7']
# c_coef_div['abs_div'] = c_coef_div['div'].abs()
c_coef_div['norm_div'] = min_max_scaler.fit_transform(c_coef_div[['div']])
#%%
standard_scaler = preprocessing.StandardScaler()
c_coef_div = mean_abs_coef_c3.merge(mean_abs_coef_c7, on = 'Aligned Position')
c_coef_div[['coef_c3','coef_c7']] = standard_scaler.fit_transform(c_coef_div[['coef_c3','coef_c7']])
c_coef_div['div'] = c_coef_div['coef_c3'] - c_coef_div['coef_c7']
# c_coef_div['abs_div'] = c_coef_div['div'].abs()
c_coef_div['norm_div'] = min_max_scaler.fit_transform(c_coef_div[['div']])

#%%

with open('2ql5.pdb', 'r') as infile,  open('c7_coef.pdb' ,'w') as outfile:
        
        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in mean_abs_coef_c7.loc[mean_abs_coef_c7['Position_c7'] == int(pdb_pos), 'WT_aa_c7']:
                    print('fuck up, aa seq doesnt align')
                
                if int(pdb_pos) not in list(mean_abs_coef_c7['Position_c7']):
                    re = str(0.00)
                else:
                    
                    re = str('%0.2f' % np.abs(mean_abs_coef_c7.loc[mean_abs_coef_c7['Position_c7'] == int(pdb_pos), 'coef_c7'])).rjust(6)
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
                if aa in mean_abs_coef_c7.loc[mean_abs_coef_c7['Position_c7'] == int(pdb_pos)-300, 'WT_aa_c7']:
                    print('fuck up, aa seq doesnt align')
                if int(pdb_pos)-300 not in list(mean_abs_coef_c7['Position_c7']):
                    re = str(0.00)
                else:
                    re = str('%0.2f' % np.abs(mean_abs_coef_c7.loc[mean_abs_coef_c7['Position_c7'] == int(pdb_pos)-300, 'coef_c7'])).rjust(6)
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

with open('2h5j.pdb', 'r') as infile,  open('c3_coef.pdb' ,'w') as outfile:
        
        
        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in mean_abs_coef_c3.loc[mean_abs_coef_c3['Position_c3'] == int(pdb_pos), 'WT_aa_c3']:
                    print('fuck up, aa seq doesnt align')
                if int(pdb_pos) not in list(mean_abs_coef_c3['Position_c3']):
                    re = str(0.00)

                else:
                    re = str('%0.2f' % np.abs(mean_abs_coef_c3.loc[mean_abs_coef_c3['Position_c3'] == int(pdb_pos), 'coef_c3'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
            
            if line[:4] == ('ATOM') and (line[21] in 'CD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                if aa in mean_abs_coef_c3.loc[mean_abs_coef_c3['Position_c3'] == int(pdb_pos), 'WT_aa_c3']:
                    print('fuck up, aa seq doesnt align')
                if int(pdb_pos) not in list(mean_abs_coef_c3['Position_c3']):
                    re = str(0.00)
                else:
                    re = str('%0.2f' % np.abs(mean_abs_coef_c3.loc[mean_abs_coef_c3['Position_c3'] == int(pdb_pos), 'coef_c3'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
           
            else: 
                outfile.write(line)
outfile.close()

#%%
with open('2ql5.pdb', 'r') as infile,  open('c7_coef_abs_div.pdb' ,'w') as outfile:
        
        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_coef_div.loc[c_coef_div['Position_c7'] == int(pdb_pos), 'WT_aa_c7']:
                    print('fuck up, aa seq doesnt align')
                
                if int(pdb_pos) not in list(c_coef_div['Position_c7']):
                    re = str(0.00)
                else:
                    re = str('%0.2f' % np.abs(c_coef_div.loc[c_coef_div['Position_c7'] == int(pdb_pos), 'norm_div'])).rjust(6)
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
                if aa in c_coef_div.loc[c_coef_div['Position_c7'] == int(pdb_pos)-300, 'WT_aa_c7']:
                    print('fuck up, aa seq doesnt align')
                if int(pdb_pos)-300 not in list(c_coef_div['Position_c7']):
                    re = str(0.00)
                else:
                    re = str('%0.2f' % np.abs(c_coef_div.loc[c_coef_div['Position_c7'] == int(pdb_pos)-300, 'norm_div'])).rjust(6)
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
with open('2h5j.pdb', 'r') as infile,  open('c3_coef_abs_div.pdb' ,'w') as outfile:
        
        
        for line in infile:

            if line[:4] == ('ATOM') and (line[21] in 'AB'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                
                if aa in c_coef_div.loc[c_coef_div['Position_c3'] == int(pdb_pos), 'WT_aa_c3']:
                    print('fuck up, aa seq doesnt align')
                if int(pdb_pos) not in list(c_coef_div['Position_c3']):
                    re = str(0.00)

                else:
                    re = str('%0.2f' % np.abs(c_coef_div.loc[c_coef_div['Position_c3'] == int(pdb_pos), 'norm_div'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
            
            if line[:4] == ('ATOM') and (line[21] in 'CD'):
                
                pdb_pos = line[23:26]              
                b = line[60:66]
                aa = three_to_one(line[17:20])
                if aa in c_coef_div.loc[c_coef_div['Position_c3'] == int(pdb_pos), 'WT_aa_c3']:
                    print('fuck up, aa seq doesnt align')
                if int(pdb_pos) not in list(c_coef_div['Position_c3']):
                    re = str(0.00)
                else:
                    re = str('%0.2f' % np.abs(c_coef_div.loc[c_coef_div['Position_c3'] == int(pdb_pos), 'norm_div'])).rjust(6)
                line = line.replace(b, re)
                outfile.write(line)
           
            else: 
                outfile.write(line)
outfile.close()

#%%
# x = c_coef_div['Aligned Position']
# y = c_coef_div['div'].rolling(3).mean()

# sns.lineplot(x = x,
#              y = y)
# plt.show()


#%%










