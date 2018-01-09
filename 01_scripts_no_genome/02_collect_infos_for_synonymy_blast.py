#!/usr/bin/env python3
"""Get sequence, locus id, position, variants

Usage:
    ./01_scripts/02_collect_infos_for_synonymy_blast.py
        wanted_loci.ids
        03_data/batch_1.sumstats.tsv
        03_data/batch_1.catalog.snps.tsv.gz
        wanted_loci.fasta
        output_file.fasta
"""

# Module
from collections import defaultdict
import gzip
import sys
import os

# Defining classes
class Fasta(object):
    """Fasta object with name and sequence
    """
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

# Defining functions
def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator
    """
    with open(input_file) as f:
        sequence = ""
        name = ""
        begun = False
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if begun:
                    yield Fasta(name, sequence)
                name = line.replace(">", "")
                sequence = ""
                begun = True
            else:
                sequence += line
        yield Fasta(name, sequence)

# Parse user input
try:
    wanted_ids_file = sys.argv[1]
    sumstats_file = sys.argv[2]
    catalog_snps_file = sys.argv[3]
    wanted_fasta_file = sys.argv[4]
    output_file = sys.argv[5]
except:
    print(__doc__)
    sys.exit(1)

# Read ids
ids = set()
for line in open(wanted_ids_file):
    l = tuple(line.strip().split())
    ids.add(l)

# Get position from sumstat
print("Get positions from sumstat")
infos = set()
for line in open(sumstats_file):
    if line.startswith("#"):
        continue

    l = line.strip().split()
    current_id = (l[1], l[3])
    if current_id in ids:
        info = (l[1], l[4])
        infos.add(info)

# Get variants from catalog.snps
print("Get variants from catalog.snps")
variants = defaultdict(list)
debug_count = 0
for line in gzip.open(catalog_snps_file, "rt"):
    if line.startswith("#"):
        continue

    l = line.strip().split()
    current_info = (l[2], l[3])

    if current_info in infos:
        variants[l[2]].append((l[3], l[6], l[7]))

# Read fasta file
print("Treat fasta file")
sequences = fasta_iterator(wanted_fasta_file)
with open(output_file, "w") as outfile:
    for s in sequences:
        locus_id = s.name.replace("locus_", "")
        infos = variants[locus_id]
        for snp in infos:
            pos = snp[0]
            var = snp[1:]

            for v in var:

                temp_seq = Fasta(s.name, s.sequence)
                temp_seq.name = temp_seq.name + "_" + pos + "_" + v
                seq = temp_seq.sequence
                temp_seq.sequence = seq[:int(pos)] + v + seq[int(pos) + 1:]
                temp_seq.write_to_file(outfile)
