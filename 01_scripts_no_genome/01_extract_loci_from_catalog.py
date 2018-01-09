#!/usr/bin/env python
"""Extract loci sequences from batch_1.catalog.tags.tsv.gz

Usage:
    ./01_scripts/01_extract_loci_from_catalog.py wanted_loci.ids catalog.tags.tsv.gz wanted_loci.fasta

wanted_loci.ids = columns 3 and 2 (in that order) from stacks 1.40 vcf output separated by a tabulation
"""

# Modules
import gzip
import sys

# Parse user input
try:
    loci_ids_file = sys.argv[1]
    catalog_file = sys.argv[2]
    loci_fasta_file = sys.argv[3]
except:
    print __doc__
    sys.exit(1)

# Main
# Read wanted ids
wanted_ids = set()
with open(loci_ids_file) as infile:
    for line in infile:
        locus_id = line.strip().split("\t")[0]
        wanted_ids.add(locus_id)

# Read catalog
print("Reading and filtering catalog")
with gzip.open(catalog_file) as cfile:
    with open(loci_fasta_file, "w") as outfile:
        for line in cfile:
            if line.startswith("#"):
                continue
            l = line.strip().split()
            current_id = l[2]
            sequence = l[8]
            if current_id in wanted_ids:
                outfile.write(">locus_" + current_id + "\n" + sequence + "\n")
