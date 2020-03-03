#!/usr/bin/python

import os
import os.path
import getopt
import sys

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["fasta_directory="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)

    input_dir = None

    for o, a in opts:
        if o == "--fasta_directory":
            input_dir = a.strip()

    assert input_dir is not None


    fasta_names = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]
    fasta_names = [f for f in fasta_names if f.endswith(".fa")]

    print "We have discovered {0} fasta files in the directory: {1}".format(str(len(fasta_names)), input_dir)
    print "Generating gtf format for each files.... Copy below output... \n\n\n\n"

    num = 1
    gtf_lines = []
    for fasta_name in fasta_names:
        seq = ""
        gene_name = None
        length = None
        with open(os.path.join(input_dir,fasta_name), "r") as fasta_file:
            fasta_f = fasta_file.readlines()
            for count, line in enumerate(fasta_f): 
                if count == 0:
                    line = line.strip().strip(">")
                    line_parts = line.split(" ")
                    gene_name = line_parts[0].strip()
                else:
                    seq += line.strip()
            #print gene_name, length, len(seq)
            chrm = "TRANS{0}".format(str(num))
            num += 1
            start = 1
            stop = len(seq)
            gtf_component = 'gene_id \"{0}\"; transcript_id \"{1}\"'.format(gene_name, gene_name+"_T_id")
            gtf_line = [chrm, "transgene_plasmid", "exon", str(start), str(stop), ".", "-", ".", gtf_component]
            print "\t".join(gtf_line)
            gtf_lines.append("\t".join(gtf_line))


    with open(os.path.join(input_dir, "all_trans.gtf"), "w") as output_file:
        for line in gtf_lines:
            output_file.write(line + "\n")

print "Complete"


'''

red-ensembl.gtf
#!genome-build GRCm38.p5
#!genome-version GRCm38
#!genome-date 2012-01
#!genome-build-accession NCBI:GCA_000001635.7
#!genebuild-last-updated 2017-01
1   ensembl_havana  gene    3205901 3671498 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2";
1   havana  transcript  3205901 3216344 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000162897"; transcript_version "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-003"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTMUST00000086625"; havana_transcript_version "1"; transcript_support_level "1";
1   havana  exon    3213609 3216344 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000162897"; transcript_version "1"; exon_number "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-003"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTMUST00000086625"; havana_transcript_version "1"; exon_id "ENSMUSE00000858910"; exon_version "1"; transcript_support_level "1";
1   havana  exon    3205901 3207317 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000162897"; transcript_version "1"; exon_number "2"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-003"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTMUST00000086625"; havana_transcript_version "1"; exon_id "ENSMUSE00000866652"; exon_version "1"; transcript_support_level "1";
1   havana  transcript  3206523 3215632 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000159265"; transcript_version "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-002"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTMUST00000086624"; havana_transcript_version "1"; transcript_support_level "1";
1   havana  exon    3213439 3215632 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000159265"; transcript_version "1"; exon_number "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-002"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTMUST00000086624"; havana_transcript_version "1"; exon_id "ENSMUSE00000863980"; exon_version "1"; transcript_support_level "1";
1   havana  exon    3206523 3207317 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000159265"; transcript_version "1"; exon_number "2"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-002"; transcript_source "havana"; transcript_biotype "processed_transcript"; havana_transcript "OTTMUST00000086624"; havana_transcript_version "1"; exon_id "ENSMUSE00000867897"; exon_version "1"; transcript_support_level "1";
1   ensembl_havana  transcript  3214482 3671498 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; tag "basic"; transcript_support_level "1";
1   ensembl_havana  exon    3670552 3671498 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; exon_number "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; exon_id "ENSMUSE00000485541"; exon_version "3"; tag "basic"; transcript_support_level "1";
1   ensembl_havana  CDS 3670552 3671348 .   -   0   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; exon_number "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; protein_id "ENSMUSP00000070648"; protein_version "4"; tag "basic"; transcript_support_level "1";
1   ensembl_havana  start_codon 3671346 3671348 .   -   0   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; exon_number "1"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; tag "basic"; transcript_support_level "1";
1   ensembl_havana  exon    3421702 3421901 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; exon_number "2"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; exon_id "ENSMUSE00000449517"; exon_version "3"; tag "basic"; transcript_support_level "1";
1   ensembl_havana  CDS 3421702 3421901 .   -   1   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; exon_number "2"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; protein_id "ENSMUSP00000070648"; protein_version "4"; tag "basic"; transcript_support_level "1";
1   ensembl_havana  exon    3214482 3216968 .   -   .   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; exon_number "3"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; exon_id "ENSMUSE00000448840"; exon_version "2"; tag "basic"; transcript_support_level "1";
1   ensembl_havana  CDS 3216025 3216968 .   -   2   gene_id "ENSMUSG00000051951"; gene_version "5"; transcript_id "ENSMUST00000070533"; transcript_version "4"; exon_number "3"; gene_name "Xkr4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000026353"; havana_gene_version "2"; transcript_name "Xkr4-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS14803"; havana_transcript "OTTMUST00000065166"; havana_transcript_version "1"; protein_id "ENSMUSP00000070648"; protein_version "4"; tag "basic"; transcript_support_level "1";
'''