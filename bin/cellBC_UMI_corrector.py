# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 22:35:38 2020

@author: heyingdong
"""
import argparse
parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
parser.add_argument('--input', required=True, help = 'Barcodes and UMI of FL reads in fastq.')
parser.add_argument('--output', required=True, help = 'Output table of barcodes.')
parser.add_argument('--tmap', required=True, help = 'Tmap file.')
parser.add_argument('-w', '--whitelist', required=True, help = 'Path of the barcode whitelist.')
parser.add_argument('--class_code', default = 'sxiypru', help = 'Code of classes specifed as unknown gene')
parser.add_argument('-t', '--threshold', required=True,type = float,help = 'Threshold of p-value for testing whether a barcode should be corrected.')
parser.add_argument('--cellranger_path', required=True, help = 'Output fasta with barcodes corrected using whitelist.')

args = parser.parse_args()

import sys
sys.path.append(args.cellranger_path+"/cellranger-cs/3.1.0/tenkit/lib/python/")
sys.path.append(args.cellranger_path+"/cellranger-cs/3.1.0/lib/python/")
import string
import array
import numpy as np
import cellranger.utils as cr_utils
import cellranger.stats as cr_stats
import tenkit.constants as tk_constants
import tenkit.seq as tk_seq
import gzip
from collections import Counter
from operator import itemgetter


def qvseq2prec(qvseq):
    return((1-10**(-(np.fromstring(qvseq, dtype=np.byte) - 33).astype(float)/10)).mean())
'''
    sumqv=0;
    for baseqv in qvseq:
        ac=1-0.1**((ord(baseqv)-33)/10)
	# qv=1-0.1**((ord(letter)-33)/10)
	sumqv+=ac
    
    qvlen=len(qvseq)
    return(sumqv/qvlen)
'''

print('Loading whitelist...')
if args.whitelist[-3:] == '.gz':
    wl = set(gzip.open(args.whitelist).read().split('\n'))
else:
    wl = set(open(args.whitelist).read().split('\n'))

print('Loading id-barcode-umi dictionary...')
id2all = {}
bc_dict_nwl = {}
wl_dict = Counter()

with open(args.input,'r') as f:
    lines = [f.readline() for i in range(4)]
    while lines[0]:
        readid = lines[0].strip()[1:]
        barcode = lines[1].strip()[:16]
        umi = lines[1].strip()[16:]
        barcodeqv = lines[3].strip()[:16]
        umiqv = lines[3].strip()[16:]
        id2all[readid] = [barcode,umi,qvseq2prec(barcodeqv),qvseq2prec(umiqv)]
        if barcode in wl:
            wl_dict.update([barcode])
            id2all[readid].append('whitelist')
            id2all[readid].append('uncorrected')
        else:
            bc_dict_nwl[readid] = [barcode,barcodeqv]
            id2all[readid].append('discarded')
            id2all[readid].append('uncorrected')
        lines = [f.readline() for i in range(4)]


print('Barcode correction...')
wl_sum = np.array(list(wl_dict.values())).sum()
wl_dict = dict(wl_dict)
for i in wl_dict:
    wl_dict[i] = float(wl_dict[i])/wl_sum

for readid in bc_dict_nwl:
    crbc = cr_stats.correct_bc_error(args.threshold,bc_dict_nwl[readid][0],bc_dict_nwl[readid][1],wl_dict)
    if crbc:
        id2all[readid][4] = 'corrected'
        id2all[readid][0] = crbc
    else:
        id2all[readid][0] = 'NA'

barcode_umi_gene2id = {}
bc_gene = {}

print('Loading gene information...')
with open(args.tmap,'r') as f:
    line = f.readline()
    line = f.readline()
    while line:
        d = line.strip().split('\t')
        if id2all[d[4]][0] != 'NA':
            gene = d[0] + '\t' + d[10]
            if d[2] in args.class_code:
                gene = gene + '\t' + d[2]           
            else:
                gene = gene + '\t' + '='
            if barcode_umi_gene2id.get((id2all[d[4]][0],id2all[d[4]][1],gene),0) == 0:
                barcode_umi_gene2id[(id2all[d[4]][0],id2all[d[4]][1],gene)] = [d[4]]
            else:
                barcode_umi_gene2id[(id2all[d[4]][0],id2all[d[4]][1],gene)].append(d[4])
            if bc_gene.get((id2all[d[4]][0],gene),0) == 0:
                bc_gene[(id2all[d[4]][0],gene)] = {}
                bc_gene[(id2all[d[4]][0],gene)][id2all[d[4]][1]] = 1
            else:
                if bc_gene[(id2all[d[4]][0],gene)].get(id2all[d[4]][1],0) == 0:
                    bc_gene[(id2all[d[4]][0],gene)][id2all[d[4]][1]] = 1
                else:
                    bc_gene[(id2all[d[4]][0],gene)][id2all[d[4]][1]] += 1
        line = f.readline()

print('UMI correction...')
correction_dict = cr_stats.correct_umis(bc_gene)
bc_gene_cr = bc_gene.copy()

cr_dict = {}
for bg, umis in correction_dict.iteritems():
    for umi, c_umi in umis.iteritems():
        bc_gene_cr[bg][c_umi] += bc_gene[bg][umi]
        bc_gene_cr[bg][umi] -= bc_gene[bg][umi]     
        cr_dict[(bg[0],umi,bg[1])] = c_umi

del bc_gene
del correction_dict

bc_umi = {}

for bg, umis in bc_gene_cr.iteritems():
    for umi, umi_count in umis.iteritems():
        if bc_umi.get((bg[0],umi),0) == 0:
            if umi_count>0:
                bc_umi[(bg[0],umi)] = {}
                bc_umi[(bg[0],umi)][bg[1]] = umi_count
        else:
            if umi_count>0:
                bc_umi[(bg[0],umi)][bg[1]] = umi_count

del bc_gene_cr

dc_dict = {}
for bu, genes in bc_umi.iteritems():
    genes_list = [[gene, count] for gene, count in genes.iteritems()]
    if len(genes_list)>1:
        genes_list = sorted(genes_list, key = itemgetter(1),reverse = True)
        if genes_list[0][1] == genes_list[1][1]:
            for i in range(len(genes_list)):
                dc_dict[(bu[0],bu[1],genes_list[i][0])] = True
        else:
            for i in range(1,len(genes_list)):
                dc_dict[(bu[0],bu[1],genes_list[i][0])] = True


for bug, readids in barcode_umi_gene2id.iteritems():
    for readid in readids:
        if dc_dict.get(bug,0) != 0:
            id2all[readid][1] = 'NA'
            id2all[readid][5] = 'discarded' 
        else:
            if cr_dict.get(bug,0) != 0:
                id2all[readid][1] = cr_dict[bug]
                id2all[readid][5] = 'corrected'

del barcode_umi_gene2id

print('Writing results...')
with open(args.output,'w') as f:
    for readid,content in id2all.iteritems():
        f.write('\t'.join([str(i) for i in [readid] + content])+'\n')
