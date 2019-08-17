#!/usr/bin/python3

import argparse
import numpy as np

from utils import out, read_gtf, read_fasta
from mapping import find_location

def run(args):
    fasta = read_fasta(args.fasta, args.logfile)
    gtf   = read_gtf(args.gtf, args.logfile)

    out("Random Indexing Test:", args.logfile)
    for i, (ID,sequence) in enumerate(fasta.items()):

        print(ID)
        if "|" in ID:
            ID = ID.split("|")[3]
        gtf_items = gtf[ID]
        if not gtf_items:
            ID = ".".join(ID.split(".")[:-1])
            gtf_items = gtf[ID]
        print(ID)

        gtf_items = gtf[ID]
        exons = [(start,end) for _,start,end,_ in gtf_items]
        print(gtf_items)
        strand = gtf_items[0][3]
        chrom  = gtf_items[0][0]

        out("Transcript ID: {}".format(ID), args.logfile)
        out("len: {:4d} sequence[:20]: {}".format(len(sequence), sequence[:20]), args.logfile)
        out("gtf_items: {}".format(gtf_items), args.logfile)

        # DEBUG: Checking that the exon lengths sum up to the length of the sequence
        total_length = 0
        for _,start,stop,_ in gtf_items:
            assert stop > start
            total_length += stop-start
        assert total_length == len(sequence), "{} != {}\n{}\n".format(total_length, len(sequence), ID, gtf_items)
        assert strand == "+" or strand == "-", "{}".format(strand)

        out("exons: {}".format(exons), args.logfile)
        for _ in range(10):
            queried_index = np.random.randint(0, len(sequence))
            genomic_index = find_location(queried_index, exons, strand, sequence)
            out("Queried Index: {}\tGenomic index: {}".format(queried_index, genomic_index), args.logfile)
            out("Transcript window around queried index: {}".format(sequence[queried_index-10:queried_index+11]), args.logfile)

        out("", args.logfile)
        if i > 20:
            exit()

    out("", args.logfile)



def main():
    parser = argparse.ArgumentParser()
    parser.set_defaults(callback=run,
                        relative_position_feature=True,
                    )
    parser.add_argument("--fasta", type=str, required=True)
    parser.add_argument("--gtf",   type=str, required=True)
    parser.add_argument("--output-log", type=str, default="log_testing.txt")
    args = parser.parse_args()
    args.logfile = open(args.output_log, "w")

    args.callback(args)


if __name__ == "__main__":
    main()
