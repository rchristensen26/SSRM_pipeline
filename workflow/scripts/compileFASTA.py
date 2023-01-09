""""
This code makes a compiled FASTA file from a list of accession numbers in batch.
Input start and end indexes to specify splice of accession number list.
"""

import os
import sys

accessionNumFile = sys.argv[1]
fasta_dir = sys.argv[2]
compiled_fasta_dir = sys.argv[3]
start_index = int(sys.argv[4])
end_index = int(sys.argv[5])

# make list of accession numbers
with open(accessionNumFile, 'r') as f:
    accession_num_file = f.read()
    accession_numbers_all = accession_num_file.split()

# splice accession numbers to create batch list
accession_list = accession_numbers_all[start_index:end_index:]

# make compiled fasta file (name includes start and end indexes)
compiled_fasta_bact_f = compiled_fasta_dir + "/" + \
                   str(start_index) + "_" + \
                   str(end_index) + "_" + \
                   "compiled_dsrAB_bact_hits.faa"

compiled_fasta_arch_f = compiled_fasta_dir + "/" + \
                   str(start_index) + "_" + \
                   str(end_index) + "_" + \
                   "compiled_dsrAB_arch_hits.faa"

# bacteria!
# check if compiled fasta already exists
if os.path.isfile(compiled_fasta_bact_f):
    print(compiled_fasta_bact_f + " already exists!")

if not os.path.isfile(compiled_fasta_bact_f):
    # iterate through accession number list
    for accession_num in accession_list:
        fasta_f_bact = fasta_dir + "/" + accession_num + "_dsrAB_bact.faa"

        if os.path.isfile(fasta_f_bact):
            os.system("cat " + fasta_f_bact + ">>" + compiled_fasta_bact_f)

    print("done with bacteria!")

# archaea!
# check if compiled fasta already exists
if os.path.isfile(compiled_fasta_arch_f):
    print(compiled_fasta_arch_f + " already exists!")

if not os.path.isfile(compiled_fasta_arch_f):
    # iterate through accession number list
    for accession_num in accession_list:
        fasta_f_arch = fasta_dir + "/" + accession_num + "_dsrAB_arch.faa"

        if os.path.isfile(fasta_f_arch):
            os.system("cat " + fasta_f_arch + ">>" + compiled_fasta_arch_f)

    print("done with archaea!")
