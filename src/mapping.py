#!/usr/bin/python3
import argparse, os
import numpy as np

from collections import defaultdict
from utils import out, read_gtf, read_bed, read_meRanCall, read_fasta

# The function for finding a mapped genome coordinate,
# Given a queried transcript index.
# Input: queried transcript index, the gtf spans, the strand, the sequence (optional)
# Output: The genomic coordinates, given in 1-based indexing.
def find_location(queried_index, coordinates, sequence=None, logfile=None):
    strand = coordinates[0][3]
    assert strand == "+" or strand == "-", "{}\n".format(strand)
    if queried_index < 0 or \
       (sequence is not None and queried_index >= len(sequence)):
        raise AssertionError("Invalid Queried Index")

    index = 0
    for coordinate in coordinates:
        _,start,stop,_ = coordinate

        difference = stop - start
        if index + difference > queried_index:
            remainder = queried_index - index
            if strand == "+":
                index = start + remainder
            else:
                # We subtract 1 here, to reconcile that we stored
                # The exon span as 0-based non-inclusive.
                # We have to subtract 1 from the end of the exon.
                index = (stop - 1) - remainder

            # We need to add 1, because GTF files have 1-based indexing.
            # And internally, we have everything in 0-based.
            return index + 1
        else:
            index = index + difference

    raise AssertionError("We should have found the queried index within the spans by now")


def run(args):
    sequences = read_fasta(args.fasta, args.logfile) if args.fasta is not None else None
    if args.gtf is not None:
        ref   = read_gtf(args.gtf, args.logfile)
    elif args.bed is not None:
        ref   = read_bed(args.bed, args.logfile)

    if ref is None:
        raise NotImplementedError

    meRanCall = read_meRanCall(args.input, args.logfile)
    outfile = open(args.output, "w")
    errfile = open(args.error_out, "w")

    header = "# <chromosome> <genomic position> <strand> <transcript ID>\n"
    outfile.write(header)

    sequence = None
    not_found_errors = 0
    not_found_error_transcripts = set()
    sequence_length_check_errors = 0
    sequence_length_check_error_transcripts = set()

    for ID, refPos, refStrand in meRanCall:
        if sequences is not None:
            sequence = sequences[ID]
        
        out("ID:{} position:{}".format(ID, refPos), args.logfile)
        if "|" in ID:
            ID = ID.split("|")[3]

        if ID not in ref:
            out("Error: Could not find ID {} in reference file".format(ID), args.logfile)
            errfile.write("Error: Could not find ID {} in reference file\n".format(ID))
            not_found_errors += 1
            not_found_error_transcripts.add(ID)
            continue
        coordinates = ref[ID]
        assert refStrand == "+"

        if sequence is not None:
            total_length = 0
            for _,start,stop,_ in coordinates:
                total_length += stop-start
            if total_length != len(sequence):
                out("Error: The lengths of the exons do not match the overall length of the sequence: exons {} != seq {}".format(total_length, len(sequence)), args.logfile)

                if len(sequence) < total_length:
                    errfile.write("Error: The lengths of the exons of ID {} do not match the overall length of the sequence: exons {} != seq {}\n".format(ID, total_length, len(sequence)))
                    sequence_length_check_errors += 1
                    sequence_length_check_error_transcripts.add(ID)
                    continue

                out("Attempting to handle Poly-A...", args.logfile)
                length_diff = len(sequence) - total_length
                last_chunk = sequence[-length_diff:]
                if last_chunk == "A"*length_diff:
                    sequence = sequence[:-length_diff]
                    out("Handled this case by cutting off Poly-A", args.logfile)
                    assert total_length == len(sequence)
                else:
                    out("Could not handle this case by cutting off Poly-A", args.logfile)
                    errfile.write("Error: The lengths of the exons of ID {} do not match the overall length of the sequence: exons {} != seq {}\n".format(ID, total_length, len(sequence)))
                    sequence_length_check_errors += 1
                    sequence_length_check_error_transcripts.add(ID)
                    continue

        strand = coordinates[0][3]
        chrom  = coordinates[0][0]
        genomic_position = find_location(refPos, coordinates)

        line = "\t".join((chrom, str(genomic_position), strand, ID)) + "\n"
        outfile.write(line)

    out("Finished writing to: {}".format(args.output), args.logfile)
    out("Number of sites whose transcript was not found in the reference file: {}".format(not_found_errors), args.logfile)
    out("Number of sites whose transcript sequence length did not equal the sum of all the exons: {}".format(sequence_length_check_errors), args.logfile)
    out("Total meRanCall lines: {}".format(len(meRanCall)), args.logfile)
    out("", args.logfile)
    out("Number of transcripts that were not found in the reference file: {}".format(len(not_found_error_transcripts)), args.logfile)
    out("Number of transcripts whose length did not equal the sum of all the exons: {}".format(len(sequence_length_check_error_transcripts)), args.logfile)
    out("Total number of transcripts read in the reference file: {}".format(len(ref)), args.logfile)

    errfile.write("Number of sites whose transcript was not found in the reference file: {}\n".format(not_found_errors))
    errfile.write("Number of sites whose transcript sequence length did not equal the sum of all the exons: {}\n".format(sequence_length_check_errors))
    errfile.write("Total meRanCall lines: {}\n".format(len(meRanCall)))
    errfile.write("\n")
    errfile.write("Number of transcripts that were not found in the reference file: {}\n".format(len(not_found_error_transcripts)))
    errfile.write("Number of transcripts whose length did not equal the sum of all the exons: {}\n".format(len(sequence_length_check_error_transcripts)))
    errfile.write("Total number of transcripts read in the reference file: {}\n".format(len(ref)))



def main():
    parser = argparse.ArgumentParser()
    parser.set_defaults(callback=run,
                        relative_position_feature=True,
                    )
    parser.add_argument("--gtf",   type=str)
    parser.add_argument("--gff",   type=str)
    parser.add_argument("--bed",   type=str)
    parser.add_argument("--fasta",   type=str)
    parser.add_argument("--input", type=str, required=True, 
                        help="Input meRanCall File is a tab-delimited file with the first 3 columns being the Transcript_ID, Position, and Strand respectively")
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    # If the user inputs a gff file instead, 
    # our read_gtf function handles either
    if args.gff is not None:
        args.gtf = args.gff

    args.output_log = os.path.splitext(os.path.basename(args.input))[0] + "_log.txt"
    args.logfile = open(args.output_log, "w") if args.output_log is not None else None
    
    args.error_out = os.path.splitext(args.output)[0] + ".err"

    args.callback(args)

if __name__ == "__main__":
    main()
