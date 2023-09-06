"""
Generates CPM data matrices for genes specified in request from NIH

Usage:

python get_gene_cpms.py \
    io: data root
    ip: path to avg vals file
    ip: path to cluster metadata file (for numcells for each cluster)
    ip: files containing requested genes list from NIH
    op: gene ppm file


Example:

python src/python/scripts/misc/for_nih/get_gene_cpms.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_avg_mtx.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/clusterSize.csv \
    /misc/for_nih/s0/mouse_GPCR_Targets.csv \
    /misc/for_nih/s1/gene_cpm_v3.csv \

Created by Mukund on 2023-08-29
"""

import csv
import sys
from produtils import dprint
import numpy as np
import re


data_root = sys.argv[1]
avg_csv_file = data_root+sys.argv[2]
clustersize_csv_file = data_root+sys.argv[3]
gene_list_file = data_root+sys.argv[4]
out_csv_file = data_root+sys.argv[5]


# get genes list
genes_list = []
# read in gene list csv file and ignore first and second rows
with open(gene_list_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        genes_list.append(row[3])

genes_list = genes_list[2:]

nClusters = 5030
nGenes = 21899

clusterSizes = []
clusterNames = []

# read clusterSize file
with open(clustersize_csv_file, 'r') as f:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(f)
    # Iterate over each row in the csv using reader object
    for idx,row in enumerate(csv_reader):
        # row variable is a list that represents a row in csv
        if (idx>0):
            clusterSizes.append(int(row[1]))
            clusterNames.append(row[0])


clusterSizes = np.array(clusterSizes)

# read avg csv file
geneNames = None
avg_mat = np.zeros(shape=(nClusters, nGenes))
with open(avg_csv_file, 'r') as f:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(f)
    # Iterate over each row in the csv using reader object
    for idx,row in enumerate(csv_reader):
        # row variable is a list that represents a row in csv
        if (idx>0):
            curcell_avg_counts = np.array([float(x) for x in row[1:]])
            dprint('avgs', idx)
            nohead_idx = idx-1
            avg_mat[nohead_idx, :] = curcell_avg_counts[:nGenes]
        if idx==0:
            geneNames = row[1:]
            geneNamesArray = np.asarray(geneNames)
            p = re.compile("(.+?)=.+$")
            geneNamesArray = np.asarray([p.search(x).group(1) for x in geneNamesArray])
            geneNames  = geneNamesArray.tolist()


# get gene indices
genesNotInGeneNames = []
geneIndices = []
chosenGeneNames = []
for idx, gene in enumerate(genes_list):
    # check if gene is in geneNames
    if gene in geneNames:
        geneIndices.append(geneNames.index(gene))
        chosenGeneNames.append(gene)
    elif gene == 'Agtr1a|Agtr1b':
        g1 = 'Agtr1a'
        g2 = 'Agtr1b'
        geneIndices.append(geneNames.index(g1))
        geneIndices.append(geneNames.index(g2))
        chosenGeneNames.append(g1)
        chosenGeneNames.append(g2)
    else:
        genesNotInGeneNames.append([idx+3,gene])



# multiply avg_mat with clusterSizes
sum_mat = avg_mat * clusterSizes[:, np.newaxis]

row_sums = sum_mat.sum(axis=1)

# # divide by row_sums
# cpm_mat = (sum_mat / row_sums[:, np.newaxis]) * 1000000

# divide each row by corresponding row_sum and multiply by 1e6
cpm_mat = np.zeros(shape=(nClusters, nGenes))
normed_counts_mat = np.zeros(shape=(nClusters, nGenes))
for idx, row in enumerate(sum_mat):
    normed_counts_mat[idx, :] = row / row_sums[idx] # get parts per unit
    cpm_mat[idx, :] = normed_counts_mat[idx, :] * 1000000 # get parts per million units

# sum across rows in avg_mat
# gene_avg = np.sum(avg_mat, axis=0)


# get ppm for each gene and append to a dataframe

# cpm_mat = avg_mat * 1000000
out_mat = np.zeros(shape=(nClusters, len(geneIndices)))

for gene_idx in geneIndices:

    # get column with index gene_idx from ppm_mat
    cur_gene_ppm = cpm_mat[:, gene_idx]
    out_mat[:, geneIndices.index(gene_idx)] = cur_gene_ppm

# take log of out_mat
# out_mat = np.log2(out_mat)

# write to csv
with open(out_csv_file, 'w') as f:
    writer = csv.writer(f)
    writer.writerow([''] + chosenGeneNames)
    for idx, row in enumerate(out_mat):
        writer.writerow([clusterNames[idx].split('=')[1]] + row.tolist())




dprint(len(genesNotInGeneNames))
dprint(genesNotInGeneNames)

dprint(np.shape(clusterSizes))
dprint(np.shape(clusterSizes[:, np.newaxis]))

for geneInfo in genesNotInGeneNames:
    print(geneInfo[0],',', geneInfo[1])

