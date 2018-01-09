#!/usr/bin/env python
"""Format synonymy_with_filter output into fasta and variant files for provean

Usage:
    ./01_scripts/01_create_provean_input.py synonymy.tsv 06_provean
"""

# Modules
import sys
import os

# Parse user input
try:
    input_file = sys.argv[1]
    output_folder = sys.argv[2]
except:
    print __doc__
    sys.exit(1)

# Functions
def find_diff(s1, s2):
    if len(s1) != len(s2):
        print "Sequence lengths are different!"
        print s1
        print s2
        return None
    
    found_diff = []
    for i in range(len(s1)):
        aa1 = s1[i]
        aa2 = s2[i]
        if aa1 != aa2:
            found_diff.append((i+1, aa1, aa2))
    
    return found_diff

# Let's go!
with open(input_file) as infile:
    for line in infile:
        if "LocusID" in line:
            continue

        l = line.strip().split("\t")

        locus = l[0]
        seq1 = l[3]
        seq2 = l[6]
        is_synonymous = l[9]

        if is_synonymous == "1":
            continue

        # Find differences
        diff_infos = find_diff(seq1, seq2)
        if not diff_infos:
            continue

        # Create files
        if diff_infos:
            infos = diff_infos[0]
            fasta_file = os.path.join(output_folder, locus + ".fasta")
            var_file = os.path.join(output_folder, locus + ".var")

            # Output sequences
            with open(fasta_file, "w") as ffile:
                ffile.write(">" + locus + "\n" + seq1 + "\n")

            # Output variant infos
            with open(var_file, "w") as vfile:
                vfile.write(infos[1] + str(infos[0]) + infos[2] + "\n")
