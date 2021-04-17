
#Parse VCF's
import numpy as d
import pandas as pd
import math
from pandas import ExcelWriter
from pandas import ExcelFile
import warnings
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import style
import statistics
from functools import reduce
import re
from pylab import *
from functools import reduce
plt.rc('font', family='serif')
plt.rc('axes', linewidth=2)
plt.rc({"fontname":"Arial"})
import allel
style.use('seaborn-muted')
import os
from glob import glob

dir_path = os.path.dirname(os.path.realpath(__file__))
filenames = glob('*Annotated.vcf')

dataframes = [allel.vcf_to_dataframe(f,fields=['calldata/GT','samples','variants/ALT','variants/CHROM','variants/FILTER', 'variants/ID', 'variants/POS', 'variants/QUAL', 'variants/REF','variants/CSQ']) for f in filenames]
print(dataframes)
dataframes_filtered = [file[file['FILTER_PASS']== True] for file in dataframes]
dataframes_filtered = [file[file['CSQ'].notnull()] for file in dataframes]
mut_type = ['missense_variant','frameshift_variant']#,'stop_lost','stop_gained']
dataframes_new = []
dataframes_symbol = []
keys = []
key=-1
for frame in dataframes_filtered:
    key = key+1
    keys.append(key)
    frame[['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','SYMBOL_SOURCE','HGNC_ID','AF','CLIN_SIG','SOMATIC','PHENO']] = frame['CSQ'].str.split('|',expand=True,)
    frame['key']=key
    frame = frame.drop_duplicates(subset='SYMBOL')
    frame = frame[['ALT_1','CHROM','POS', 'QUAL', 'REF','Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','SYMBOL_SOURCE','HGNC_ID','AF','CLIN_SIG','SOMATIC','PHENO','key']].copy()
    frame = frame[frame['Consequence'].isin(mut_type)]
    dataframes_new.append(frame)
    dataframes_symbol.append(frame['SYMBOL'])
dataframes_new = pd.concat(dataframes_new)
names_list = reduce(lambda x, y: pd.merge(x, y, on = 'SYMBOL'), dataframes_symbol)
names = list(names_list['SYMBOL'])
final_frame  = dataframes_new[dataframes_new['SYMBOL'].isin(names)]
final_frame.to_csv('All_dbSNP_4-9-21.csv')



