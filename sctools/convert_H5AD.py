#!/usr/bin/python3

## Convert HD5A to cellranger style output
## Usage: python3 convert_H5AD.py --input h5ad_file.h5ad
## Author: @neil.alistair.robertson@hotmail.co.uk


import scanpy as sc
import numpy as np
import scipy.sparse as sparse
import scipy.io as scio
import gzip
import getopt, sys, os

OUTPUT_FILE_SUFFIXES =(".mtx", ".barcodes.tsv", ".features.tsv")
options, remainder = getopt.getopt(sys.argv[1:], 'i', ['input='])

input_filename = None

for opt, arg in options:
    if opt in ("-i", "--input"):
        input_filename = arg.strip()

assert input_filename != None

print("-> Reading HD5A file: {0}".format(input_filename))
scanpy_object = sc.read(input_filename)

obs = scanpy_object.obs.to_dict()
var = {"Gene" : np.array(scanpy_object.var.index)}

print("-> Is file in memory bonafide sparse array?: {0}".format(str(sparse.issparse(scanpy_object.X))))

print("-> Writing file to mtx format...")
scio.mmwrite(target = input_filename, a = scanpy_object.X)

print("-> Now extracting and writing metadata...")
scanpy_object.obs.to_csv(input_filename + ".barcodes.tsv", header = True, sep = "\t")

with open(input_filename + ".features.tsv", "w") as features_out:
    for feature in np.array(scanpy_object.var.index):
        features_out.write(feature.strip() + "\n")

dim_reduction = scanpy_object.uns

print("-> Compressing files...")
for suffix in OUTPUT_FILE_SUFFIXES:
    with open(input_filename + suffix, "r") as f_in:
        with gzip.open(input_filename + suffix + ".gz", 'wt') as f_out:
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()

print("-> Removing uncompressed files...")
for suffix in OUTPUT_FILE_SUFFIXES:
    os.remove(input_filename + suffix)

print("-> Complete! Attempt to load with scater (code below)...\n\n")
print("#################")
print("##### RCODE #####")
print("#################")

## END ##

print('''require("Matrix")
require("SingleCellExperiment")
require("scater")

barcodes <- read.table("{0}.barcodes.tsv.gz", sep = "\t", header = T)
features <- read.table("{1}.features.tsv.gz", sep = "\t")
mtx <- Matrix::readMM("{2}.mtx.gz") ###load up matrix.mtx
mtx <- t(as.matrix(mtx))

rownames(mtx) <- make.unique(as.character(features[,1])) ####if second column is genenames then we want this one, and also uniquing any duplicates, which is essentially adding .1 or .2 to duplicates 
colnames(mtx) <- make.unique(as.character(barcodes$index))

anno <- data.frame(barcodes)

SCE <- SingleCellExperiment(assays = list(counts = mtx), colData = anno)
rm(list=ls()[!(ls() %in% c("SCE"))])
SCE

logcounts(SCE) <- log2(as.matrix(counts(SCE))+1)
normalize(SCE, log = TRUE)
cpm(SCE) <- log2(calculateCPM(SCE) + 1)

png("{3}.tsne.png")
plotTSNE(SCE)
dev.off()

saveRDS(SCE, "{4}.SingleCellExperiment.Rds")
'''.format(input_filename,input_filename,input_filename,input_filename,input_filename))


