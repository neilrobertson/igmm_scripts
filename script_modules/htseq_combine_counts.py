#!/usr/bin/python

import os
import getopt
import sys

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["input=","output="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)

    for o, a in opts:
        if o == "--input":
            input_dir = a.strip()
        if o == "--output":
            output_dir = a.strip()

    assert input_dir is not None
    assert output_dir is not None

    genes_dict = {}
    samples_counts_dict_list = []

    sample_names = [dI for dI in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir,dI))]

    for sample in sample_names:
        sample_filename = os.path.join(input_dir, sample, "accepted_hits.htseq-counts.txt")
        with open(sample_filename, "r") as sample_file:
            sample_count_dict = {}
            sample_tpm_dict = {}
            for line in sample_file:
                gene, counts = line.strip().split("\t")
                genes_dict[gene] = True

                sample_count_dict[gene] = counts

            samples_counts_dict_list.append(sample_count_dict)

    counts_out_file = open(os.path.join(output_dir, "combined.htseq-counts.tsv"),"w")

    # does the output
    counts_out_file.write("ID" + "\t" + "\t".join(sample_names) + "\n")

    for gene in genes_dict:
        counts_out_list = [gene]
        for sample_count_dict in samples_counts_dict_list:
            if gene in sample_count_dict:
                counts_out_list.append(sample_count_dict[gene])
            else:
                counts_out_list.append("0")
        counts_out_file.write("\t".join(counts_out_list) + "\n")

