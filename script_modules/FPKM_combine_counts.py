#!/usr/bin/python

import os
import getopt
import sys

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["FPKM_directory="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)

    input_dir = None

    for o, a in opts:
        if o == "--FPKM_directory":
            input_dir = a.strip()

    assert input_dir is not None

    genes_dict = {}
    samples_counts_dict_list = []

    sample_names = [dI for dI in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir,dI))]

    for sample in sample_names:
        sample_filename = os.path.join(input_dir, sample, "genes.fpkm_tracking")
        with open(sample_filename, "r") as sample_file:
            sample_count_dict = {}
            for count, line in enumerate(sample_file):
                if count == 0:
                    pass
                else:
                    tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, locus, length, coverage, FPKM, FPKM_conf_lo, FPKM_conf_hi, FPKM_status = line.strip().split("\t")
                    genes_dict[gene_id] = gene_short_name
                    sample_count_dict[gene_id] = FPKM

            samples_counts_dict_list.append(sample_count_dict)

    counts_out_file = open(os.path.join(input_dir, "FPKM_matrix.csv"),"w")

    # does the output
    counts_out_file.write("EnsemblID\tGene_Name\t" + "\t".join(sample_names) + "\n")

    for gene in genes_dict:
        out_list = [gene]
        out_list.append(genes_dict[gene])
        for sample_count_dict in samples_counts_dict_list:
            if gene in sample_count_dict:
                out_list.append(sample_count_dict[gene])
            else:
                out_list.append("0")
        counts_out_file.write("\t".join(out_list) + "\n")


