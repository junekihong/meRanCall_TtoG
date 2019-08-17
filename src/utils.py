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
    for line in open(filename, "r"):
        if line[0] == "#":
            continue
        line = line.strip()
        
        transcript_ID_line = line.split(";")        
        line = line.split()

        # Take the exons only. Skip over the transcript lines.
        # We are going to walk across the exons.
        if line[2] != "exon":
            continue

        transcript_ID = [x for x in transcript_ID_line if "transcript_id" in x]
        if not transcript_ID:
            continue
        transcript_ID = transcript_ID[-1][len("transcript_id")+1:].strip().strip("\"")


        # The GTF file contains indices that are 1-based, inclusive
        # When we read them in, we read them in as 0-based, exclusive
        chrom         = line[0]
        #transcript_ID = line[11].strip(";").strip("\"")
        start, stop   = int(line[3])-1, int(line[4])
        strand        = line[6]
        assert strand == "+" or strand == "-"

        item = (chrom,start,stop,strand)
        gtf[transcript_ID].append(item)

    num_transcripts      = len(gtf)
    num_genome_locations = sum([len(x) for x in gtf.values()])
    out("...\tRead {} Transcripts and {} Genome Locations".format(num_transcripts, num_genome_locations), logfile)

    for transcript_ID in gtf:
        strand = gtf[transcript_ID][0][3]
        gtf[transcript_ID] = sorted(list(gtf[transcript_ID]), key=lambda x:x[1], reverse=strand=="-")
    return gtf

def read_bed(filename, logfile=None):
    out("Reading Bed file: {}".format(filename), logfile, end="")
    bed = defaultdict(list)

    for line in open(filename, "r"):
        if line[0] == "#":
            continue
        line = line.strip().split()

        # The Bed file contains indices that are 0-based, exclusive
        # When we read them in, we read them in as 0-based, exclusive
        chrom         = line[0]
        transcript_ID = line[3]
        start, stop   = int(line[1]), int(line[2])
        strand        = line[5]
        blockSizes    = tuple([int(x) for x in line[10].split(",")[:-1]])
        blockStarts   = tuple([int(x) for x in line[11].split(",")[:-1]])
        for bStart,bSize in zip(blockStarts, blockSizes):
            item = (chrom, start+bStart, start+bStart+bSize, strand)
            bed[transcript_ID].append(item)

    num_transcripts      = len(bed)
    num_genome_locations = sum([len(x) for x in bed.values()])
    out("...\tRead {} Transcripts and {} Genome Locations".format(num_transcripts, num_genome_locations), logfile)

    for transcript_ID in bed:
        strand = bed[transcript_ID][0][3]
        bed[transcript_ID] = sorted(list(bed[transcript_ID]), key=lambda x:x[1], reverse=strand=="-")
    return bed



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
    out("...\tRead {} Lines".format(len(result)), logfile)
    return result
