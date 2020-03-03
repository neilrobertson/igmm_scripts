#!/usr/bin/python

import getopt
import sys
import os

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["input=","output="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)


    path_list = []
    out_file = None
    chromosomes = {}

    sample = []
    input_reads = []
    uniquely_mapped_reads = []
    uniquely_mapped_read_rate = []
    number_of_splices = []
    mismatch_rate_per_base = []
    deletion_rate_per_base = []
    insertion_rate_per_base = []
    multi_mapping_reads = []
    multi_mapping_reads_rate = []
    too_many_loci_reads = []
    too_many_loci_reads_rate = []
    unmapped_too_many_mismatches_rate = []
    unmapped_too_short_rate = []
    unmapped_other_rate = []

    for o, a in opts:
        if o == "--input":
            path_root = a.strip()
        if o == "--output":
            out_filename = a.strip()

    assert path_root
    assert out_filename

    alignment_directories = [dI for dI in os.listdir(path_root) if os.path.isdir(os.path.join(path_root,dI))]

    for alignment_replicate in alignment_directories:

        file_path = "/".join([path_root, alignment_replicate, "Log.final.out"])
        print "Gathering statistics for file {0}".format(alignment_replicate)
        with open(file_path, "r") as alignment_file_summary:
            alignment_file = alignment_file_summary.readlines()

            sample.append(alignment_replicate)
            input_reads.append(alignment_file[5].rstrip().split("\t")[1])
            uniquely_mapped_reads.append(alignment_file[8].rstrip().split("\t")[1])
            uniquely_mapped_read_rate.append(alignment_file[9].rstrip().split("\t")[1])
            multi_mapping_reads.append(alignment_file[23].rstrip().split("\t")[1])
            multi_mapping_reads_rate.append(alignment_file[24].rstrip().split("\t")[1])
            too_many_loci_reads.append(alignment_file[25].rstrip().split("\t")[1])
            too_many_loci_reads_rate.append(alignment_file[26].rstrip().split("\t")[1])
    	    unmapped_too_many_mismatches_rate.append(alignment_file[28].rstrip().split("\t")[1])
    	    unmapped_too_short_rate.append(alignment_file[29].rstrip().split("\t")[1])
    	    unmapped_other_rate.append(alignment_file[30].rstrip().split("\t")[1])

            number_of_splices.append(alignment_file[11].rstrip().split("\t")[1])
            mismatch_rate_per_base.append(alignment_file[17].rstrip().split("\t")[1])
            deletion_rate_per_base.append(alignment_file[19].rstrip().split("\t")[1])
            insertion_rate_per_base.append(alignment_file[21].rstrip().split("\t")[1])


    print "Writing statistics output..."
    with open(os.path.join(out_filename, "alignment_stats.csv"), "w") as out_file:
        out_file.write("\t".join(["sample", "input_reads", "uniquely_mapped_reads", "uniquely_mapped_reads_rate","multi_mapping_reads","multi_mapping_reads_rate","reads_mapped to_too_many_loci","reads_mapped to_too_many_loci_rate","unmapped_reads_too_many_mismatches_rate","unmapped_reads_too_short_rate","unmapped_reads_other_rate", "number_of_splices","mismatch_rate_per_base", "deletion_rate_per_base", "insertion_rate_per_base"]) + "\n")
        for index in range(0,len(sample)):
            out_file.write("\t".join([sample[index],input_reads[index],uniquely_mapped_reads[index],uniquely_mapped_read_rate[index],multi_mapping_reads[index],multi_mapping_reads_rate[index],too_many_loci_reads[index],too_many_loci_reads_rate[index],unmapped_too_many_mismatches_rate[index],unmapped_too_short_rate[index],unmapped_other_rate[index],number_of_splices[index],mismatch_rate_per_base[index],deletion_rate_per_base[index],insertion_rate_per_base[index]]) + "\n")
