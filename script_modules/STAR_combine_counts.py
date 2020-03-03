#!/usr/bin/python

import os
import getopt
import sys

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["input=","output=","strand="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)

    for o, a in opts:
        if o == "--input":
            input_dir = a.strip()
        if o == "--output":
            output_dir = a.strip()
        if o == "--strand":
            strand = a.strip()

    assert input_dir is not None
    assert output_dir is not None
    assert strand in ("0","1","2")

    genes_dict = {}
    genes_dict["N_total"] = True
    samples_counts_dict_list = []

    sample_names = [dI for dI in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir,dI))]

    for sample in sample_names:
        sample_filename = os.path.join(input_dir, sample, "ReadsPerGene.out.tab")
        with open(sample_filename, "r") as sample_file:
            sample_count_dict = {}
            sample_total = 0
            for line in sample_file:
                gene, all_counts, stranded, reverse_stranded = line.strip().split("\t")
                genes_dict[gene] = True

                if strand == "0":
                    sample_count_dict[gene] = all_counts
                    sample_total += int(all_counts  )
                if strand == "1":
                    sample_count_dict[gene] = stranded
                    sample_total += int(stranded)
                if strand == "2":
                    sample_count_dict[gene] = reverse_stranded
                    sample_total += int(reverse_stranded)

            sample_count_dict["N_total"] = str(sample_total)
            samples_counts_dict_list.append(sample_count_dict)

    counts_out_file = open(os.path.join(output_dir, "read_counts.csv"),"w")
    summary_out_file = open(os.path.join(output_dir, "read_counts_summary.csv"),"w")

    # does the output
    counts_out_file.write("\t" + "\t".join(sample_names) + "\n")
    summary_out_file.write("\t" + "\t".join(sample_names) + "\n")

    for gene in genes_dict:

        out_list = [gene]
        for sample_count_dict in samples_counts_dict_list:
            if gene in sample_count_dict:
                out_list.append(sample_count_dict[gene])
            else:
                out_list.append("0")

        if gene == "N_unmapped" or gene == "N_noFeature" or gene == "N_multimapping" or gene == "N_ambiguous" or gene == "N_total":
            summary_out_file.write("\t".join(out_list) + "\n")
        else:
            counts_out_file.write("\t".join(out_list) + "\n")


