# meRanCall_TtoG


Please unzip the .tar.gz files before running. 
This will provide a GTF file and a FASTA file, used for the example programs.

Usage: 
python3 src/mapping.py --gtf merged_transcripts.gtf --input meRanCall/sample_meRanCall.txt --output output/sample_meRanCall.genomic-coord.txt

You can also see example_run.sh as well, for an example usage.

You can run our tests by running
python3 src/mapping_test.py --fasta merged_transcripts.fa --gtf merged_transcripts.gtf