#!/usr/bin/env bash
set -euo pipefail

### MMSeqs2 Annotation ###
# This script uses MMSeqs2 to attribute UniRef100 ids to the sequences in the fasta file

# Edit this to your project directory
PROJ_DIR="/mnt/ssd/elliot"

# Create project directory
mdkir "$PROJ_DIR"
mdkir "$PROJ_DIR"/tmp


# First, download the UniRef100 database
mmseqs databases UniRef100 "$PROJ_DIR"/uniref100 "$PROJ_DIR"/tmp

# Then, create sequence database for query and target
mmseqs createdb "$PROJ_DIR"/uniref100 "$PROJ_DIR"/uniref100DB
mmseqs createdb "$PROJ_DIR"/panproteome-cacnes.fasta "$PROJ_DIR"/panproteome-cacnesDB

# Then, create index for target database
mmseqs createindex "$PROJ_DIR"/panproteome-cacnesDB "$PROJ_DIR"/tmp

# Then, search the query database against the target database
    # mmseqs search queryDB targetDB resultDB tmp
    # -s is the sensitivity
mmseqs search --threads 32 -s 6.0 --disk-space-limit 750G \
    "$PROJ_DIR"/panproteome-cacnesDB \
    "$PROJ_DIR"/uniref100DB \
    "$PROJ_DIR"/resultDB \
    "$PROJ_DIR"/tmp \
    

# Then, convert the output to BLAST tsv format
    # mmseqs convertalis queryDB targetDB resultDB resultDB.m8
mmseqs convertalis --threads 32 \
    "$PROJ_DIR"/panproteome-cacnesDB \
    "$PROJ_DIR"/uniref100DB \
    "$PROJ_DIR"/resultDB \
    "$PROJ_DIR"/resultDB.m8 \

# parametre pour search
# mmseqs search --help
    # --taxon-list STR
    # Taxonomy ID, possibly multiple values separated by ',' []

# https://github.com/soedinglab/MMseqs2/releases/tag/14-7e284
    # --taxon-list parameter understands expressions.
    # E.g. get all bacterial and human sequences --taxon-list "2||9606"