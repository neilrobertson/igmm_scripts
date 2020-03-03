#!/usr/bin/python

#############################################################
#### get_10xATAC_cluster_bam.py #############################
#### Author: Neil Robertson #################################
#### Generates cluster specific bams from cohort level ######
#############################################################



## Requires tqdm to make progress bar
import os
import sys
import getopt
import pysam
import csv
from tqdm import tqdm



# this is the graph based clustering file: './aggr/outs/analysis/clustering/graphclust/clusters.csv'
clusters_file = None
# This is the sample sheet used to aggregate single replicates in cellranger-atac: 'sample_sheet.csv'
sample_sheet = None

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["clusters_file=","sample_sheet="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)

for o, a in opts:
    if o == "--clusters_file":
        clusters_file = a.strip()
    if o == "--sample_sheet":
        sample_sheet = a.strip()

assert clusters_file != None
assert sample_sheet != None



cluster_dict = {}
with open(clusters_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    #skip header
    header = next(csv_reader)
    for row in csv_reader:
        cluster_dict[row[0]] = row[1]

clusters = set(x for x in cluster_dict.values())
reps = set(x.strip().split("-")[-1] for x in cluster_dict.keys())

rep_cluster_dict = {}
for i in reps:
    rep_list = [x.strip().split('-')[0] for x in cluster_dict.keys() if x.strip().split('-')[-1] == i]
    tmp_dict = {}
    for j in rep_list:
        tmp_dict["{0}-1".format(j)] = "{0}-{1}".format(j, i)
    rep_cluster_dict[i] = tmp_dict

fin_dict = {}
name_dict = {}
with open(sample_sheet) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    #skip header
    header = next(csv_reader)
    for i, row in enumerate(csv_reader):
        rep_name = row[0].strip()
        rep_bam = row[1].strip().replace("fragments.tsv.gz", "possorted_bam.bam") 
        fin = pysam.AlignmentFile(rep_bam, "rb")
        fin_dict[str(i + 1)] = fin
        name_dict[str(i + 1)] = rep_name

# open the number of bam files as the same number of clusters, and map the out file handler to the cluster id, write to a bam with wb
fouts_dict = {}
for cluster in clusters:
    fout = pysam.AlignmentFile("cluster" + cluster + ".bam", "wb", template = fin)
    fouts_dict[cluster] = fout


for i in reps:
    print "Working on replicate {0}: {1}".format(i, name_dict[i])
    fin = fin_dict[i]
    for read in tqdm(fin):
        tags = read.tags
        CB_list = [ x for x in tags if x[0] == "CB"]
        if CB_list:
            cell_barcode = CB_list[0][1]
        # the bam files may contain reads not in the final clustered barcodes
        # will be None if the barcode is not in the clusters.csv file
        else: 
            continue
        rep_barcode = rep_cluster_dict[i].get(cell_barcode)
        cluster_id = cluster_dict.get(rep_barcode)
        if cluster_id:
            fouts_dict[cluster_id].write(read)

## do not forget to close the files
for fin in fin_dict.values():
    fin.close()

for fout in fouts_dict.values():
    fout.close()

