
#Parse VCF's
import numpy as d
import pandas as pd
import math
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy
from scipy import stats
from scipy.stats import pearsonr
import warnings
from scipy.stats import ttest_ind
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import style
import statistics
import seaborn as sns
from functools import reduce
from scipy import stats
import re
from adjustText import adjust_text
from pylab import *
from functools import reduce
plt.rc('font', family='serif')
plt.rc('axes', linewidth=2)
plt.rc({"fontname":"Arial"})
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')

import allel
style.use('seaborn-muted')
import os

from glob import glob
dir_path = os.path.dirname(os.path.realpath(__file__))
filenames = glob('*_fields.vcf')

dataframes = [allel.vcf_to_dataframe(f,fields=['calldata/GT','samples','variants/ALT','variants/CHROM','variants/FILTER', 'variants/ID', 'variants/POS', 'variants/QUAL', 'variants/REF','variants/CSQ']) for f in filenames]
print(dataframes)
#Treated = Treated_frame[(Treated_frame.FILTER_PASS == True)]
dataframes_filtered = [file[file['FILTER_PASS']== True] for file in dataframes]
dataframes_filtered = [file[file['CSQ'].notnull()] for file in dataframes]
mut_type = ['missense_variant','frameshift_variant']#,'stop_lost','stop_gained']
#print(dataframes_filtered[0].CSQ)
#Get the gene names and ensembl genes, aggergate to a list
#dataframe_new = pd.DataFrame(columns = ['ALT_1','CHROM','POS', 'QUAL', 'REF','Gene','SYMBOL','Consequence','Protein_position','Existing_variation','Amino_acids','AF','CLIN_SIG','PHENO'])
dataframes_new = []
dataframes_symbol = []
keys = []
key=-1
for frame in dataframes_filtered:
    key = key+1
    keys.append(key)
    #new_frame = pd.DataFrame(columns = ['Gene','SYMBOL','Consequence','Protein_position','Existing_variation','Amino_acids','AF','CLIN_SIG','PHENO'])
    frame[['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','SYMBOL_SOURCE','HGNC_ID','AF','CLIN_SIG','SOMATIC','PHENO']] = frame['CSQ'].str.split('|',expand=True,)
    frame['key']=key
    #new_frame[['Ensembl','Name','Consequence','PositionProtein','Co-locatedVariant','AF','CLIN_SIG','PHENO']] = frame['CSQ'].str.split('|',n=8,expand=True)
    frame = frame.drop_duplicates(subset='SYMBOL')
    frame = frame[['ALT_1','CHROM','POS', 'QUAL', 'REF','Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','SYMBOL_SOURCE','HGNC_ID','AF','CLIN_SIG','SOMATIC','PHENO','key']].copy()
    frame = frame[frame['Consequence'].isin(mut_type)]
    #print(frame['SYMBOL'])
    #new_frame.index = new_frame['Ensembl']
    dataframes_new.append(frame)
    dataframes_symbol.append(frame['SYMBOL'])
#print(len(dataframes_symbol))
dataframes_new = pd.concat(dataframes_new)
#names_list = pd.merge([x['SYMBOL'] for x in dataframes_filtered])
#print(names_list)
#names_list = names_list.drop_duplicates()
#names_list = pd.concat([x.set_index('SYMBOL') for x in dataframes_symbol], axis=1, join='inner')
names_list = reduce(lambda x, y: pd.merge(x, y, on = 'SYMBOL'), dataframes_symbol)
names = list(names_list['SYMBOL'])
print(names)
print(dataframes_new.head())
#final_frame = pd.concat([x.set_index('SYMBOL') for x in dataframes_filtered], axis=1, join='inner')
#names_list = pd.concat([x.set_index('SYMBOL') for x in dataframes_filtered], axis=1, join='inner')
#final_frame = reduce(lambda x, y: pd.merge(x, y, on = 'SYMBOL'), dataframes_filtered)
#names_list = pd.concat(dataframes_new,axis=0, join='inner')
final_frame  = dataframes_new[dataframes_new['SYMBOL'].isin(names)]
print(len(final_frame))
#print(names_list)
final_frame.to_csv('All_dbSNP_4-9-21.csv')



