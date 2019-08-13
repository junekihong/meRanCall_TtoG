#!/usr/bin/python3

import sys
from collections import defaultdict
from Bio import SeqIO

def out(string, logfile, end="\n"):
    print(string, end=end)
    sys.stdout.flush()
    if logfile is not None:
        logfile.write("{}{}".format(string, end))
        logfile.flush()

def read_gtf(filename, logfile=None):
    out("Reading gtf file: {}".format(filename), logfile, end="")
    gtf = defaultdict(list)
    debug_dict = {}

    for line in open(filename, "r"):
        if line[0] == "#":
            continue
        line = line.strip().split()

        # Take the exons only. Skip over the transcript lines.
        # We are going to walk across the exons.
        if line[2] != "exon":
            continue

        # The GTF file contains indices that are 1-based, inclusive
        # When we read them in, we read them in as 0-based, exclusive
        chrom         = line[0]
        transcript_ID = line[11].strip(";").strip("\"")
        start, stop   = int(line[3])-1, int(line[4])
        strand        = line[6]
        item = (chrom,start,stop,strand)

        gtf[transcript_ID].append(item)

    num_transcripts      = len(gtf)
    num_genome_locations = sum([len(x) for x in gtf.values()])
    out("...\tRead {} Transcripts and {} Genome Locations".format(num_transcripts, num_genome_locations), logfile)

    for transcript_ID in gtf:
        gtf[transcript_ID] = sorted(list(gtf[transcript_ID]))
    return gtf

def read_fasta(filename, logfile=None):
    out("Reading Fasta File: {}".format(filename), logfile, end="")
    items = SeqIO.parse(filename, "fasta")
    fasta = {}
    for item in items:
        fasta[item.id] = item.seq
    out("...\tRead {} Transcripts".format(len(fasta)), logfile)
    return fasta

def read_meRanCall(filename, logfile=None):
    out("Reading meRanCall File: {}".format(filename), logfile, end="")

    result = []
    for line in open(filename, "r"):
        if line[0] == "#":
            continue
        line = line.strip().split()

        # The meRanCall file uses 1-based indexing.
        # We read it in as 0-based indexing.
        ID        = line[0]
        refPos    = int(line[1]) - 1
        refStrand = line[2]

        item = (ID, refPos, refStrand)
        result.append(item)
    return result
