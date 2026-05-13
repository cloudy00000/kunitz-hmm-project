#!/bin/bash

# Kunitz HMM validation pipeline
# This script assumes that the following files already exist:
# - fasta.seq
# - kunitz.hmm alignment input or alignment file for hmmbuild
# - ./Uniprot/Positive.fasta
# - ./Uniprot/Negative.fasta
# - performance.py


# 1. We go to project directory
cd ~/bioinf/bioinf_project/

# Build the Kunitz HMM from the alignment
hmmbuild kunitz.hmm fasta.seq

# We search the positive and negative UniProt datasets with hmmsearch
hmmsearch --max --noali --tblout positive_kunitz.tbl -Z 1000 kunitz.hmm ./Uniprot/Positive.fasta
hmmsearch --max --noali --tblout negative_kunitz.tbl -Z 1000 kunitz.hmm ./Uniprot/Negative.fasta

# We count the reported HMMER hits
grep -v '^#' positive_kunitz.tbl | awk '{print $1,$8}' | wc
grep -v '^#' negative_kunitz.tbl | awk '{print $1,$8}' | wc

# We create a labelled match files
# Format: protein_id    evalue    true_label
# true_label = 1 for known Kunitz positives
# true_label = 0 for known negatives
grep -v '^#' positive_kunitz.tbl | awk '{print $1"\t"$8"\t1"}' > positive_kunitz.match
grep -v '^#' negative_kunitz.tbl | awk '{print $1"\t"$8"\t0"}' > negative_kunitz.match

# We extract all negative IDs from the original negative FASTA file
grep '^>' ./Uniprot/Negative.fasta | awk '{print $1}' | tr -d '>' | sort > negative_kunitz.ids

# We extract negative IDs that had HMMER hits
awk '{print $1}' negative_kunitz.match | sort > negative_kunitz_match.ids

# We create file of negative proteins with no HMMER hit
# These are assigned an artificial bad E-value of 100 and true label 0
comm -23 <(sort negative_kunitz.ids) <(sort negative_kunitz_match.ids) \
| awk '{print $1"\t100\t0"}' > negative_kunitz.nonmatch

# We combine matched and nonmatched negatives, then shuffle
cat negative_kunitz.match negative_kunitz.nonmatch | sort -R > negative_kunitz.tot.match

# We check dataset counts
echo "Original FASTA counts:"
grep -c '^>' ./Uniprot/Positive.fasta
grep -c '^>' ./Uniprot/Negative.fasta

echo "Evaluation file counts:"
wc -l positive_kunitz.match
wc -l negative_kunitz.match
wc -l negative_kunitz.nonmatch
wc -l negative_kunitz.tot.match

# Positives are split into two validation sets
head -n 199 positive_kunitz.match > kunitz_set_1.txt
tail -n 199 positive_kunitz.match > kunitz_set_2.txt

# Negatives are added to the two validation sets
head -n 287115 negative_kunitz.tot.match >> kunitz_set_1.txt
tail -n 287115 negative_kunitz.tot.match >> kunitz_set_2.txt

# We check the validation set sizes and label distributions
echo "Validation set sizes:"
wc -l kunitz_set_1.txt
wc -l kunitz_set_2.txt

echo "Validation set 1 labels:"
awk '{print $3}' kunitz_set_1.txt | sort | uniq -c

echo "Validation set 2 labels:"
awk '{print $3}' kunitz_set_2.txt | sort | uniq -c

# We run performance evaluation across E-value thresholds
for i in $(seq 1 15); do
    python3 performance.py kunitz_set_1.txt 1e-$i
done > kunitz_set_1.results

for i in $(seq 1 15); do
    python3 performance.py kunitz_set_2.txt 1e-$i
done > kunitz_set_2.results

# We preform a final inspection of the results
echo "Performance results for validation set 1:"
cat kunitz_set_1.results

echo "Performance results for validation set 2:"
cat kunitz_set_2.results

# We make it execyutable and run it
#chmod +x kunitz_pipeline.sh
#./kunitz_pipeline.sh