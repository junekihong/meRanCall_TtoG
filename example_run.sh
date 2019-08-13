#!/bin/sh

time python3 src/mapping.py --gtf merged_transcripts.gtf --input meRanCall/sample_meRanCall.txt --output output/sample_meRanCall.genomic-coord.txt
