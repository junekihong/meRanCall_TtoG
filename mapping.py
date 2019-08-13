#!/usr/bin/python3

import argparse
import numpy as np

from utils import out, read_gtf, read_meRanCall

# The function for finding a mapped genome coordinate,
# Given a queried transcript index.
# Input: queried transcript index, the gtf spans, the strand, the sequence (optional)
# Output: The genomic coordinates, given in 1-based indexing.
def find_location(queried_index, spans, strand, sequence=None, logfile=None):
    assert strand == "+" or strand == "-"
    if queried_index < 0 or \
       (sequence is not None and queried_index >= len(sequence)):
        raise AssertionError("Invalid Queried Index")

    if strand == "-":
        spans = reversed(spans)

    index = 0
    for span in spans:
        difference = span[1] - span[0]
        if index + difference > queried_index:
            remainder = queried_index - index
            if strand == "+":
                coordinate = span[0] + remainder
            else:
                # We subtract 1 here, to reconcile that we stored
                # The exon span as 0-based non-inclusive.
                # We have to subtract 1 from the end of the exon.
                coordinate = (span[1] - 1) - remainder

            # We need to add 1, because GTF files have 1-based indexing.
            # And internally, we have everything in 0-based.
            return coordinate + 1
        else:
            index = index + difference
    raise AssertionError("We should have found the queried index within the spans by now")

def run(args):
    gtf       = read_gtf(args.gtf, args.logfile)
    meRanCall = read_meRanCall(args.input, args.logfile)
    outfile = open(args.output, "w")
    #header = "<chromosome> <genomic position> <strand>\n"
    #outfile.write(header)

    for ID, refPos, refStrand in meRanCall:
        coordinates = gtf[ID]
        assert refStrand == "+"
        
        exons  = [(start,end) for _,start,end,_ in coordinates]
        strand = coordinates[0][3]
        chrom  = coordinates[0][0]
        genomic_position = find_location(refPos, exons, strand)

        line = "\t".join((chrom, str(genomic_position), strand)) + "\n"
        outfile.write(line)
    out("Finished writing to: {}".format(args.output), args.logfile)

def main():
    parser = argparse.ArgumentParser()
    parser.set_defaults(callback=run,
                        relative_position_feature=True,
                    )
    parser.add_argument("--gtf",   type=str, required=True)
    parser.add_argument("--input", type=str, required=True, 
                        help="Input meRanCall File is a tab-delimited file with the first 3 columns being the Transcript_ID, Position, and Strand respectively")
    parser.add_argument("--output", type=str, required=True)

    #parser.add_argument("--output-log", type=str, default="log.txt")
    parser.add_argument("--output-log", type=str, default=None)

    args = parser.parse_args()

    args.logfile = open(args.output_log, "w") if args.output_log is not None else None
    args.callback(args)

if __name__ == "__main__":
    main()
