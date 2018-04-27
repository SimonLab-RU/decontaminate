from __future__ import print_function, division

import argparse
import os
import re
from collections import defaultdict, Counter
from itertools import combinations
import numpy as np
import Levenshtein

def fastq_iter(path):
    """ Quickly iterate though a fastq file (avoids bioseq dependency).
    """
    with open(path, 'r') as f:
        while True:
            try:
                header = next(f).rstrip()
                seq = next(f).rstrip()
                _ = next(f)
                qual = next(f).rstrip()
                yield header, seq, qual
            except StopIteration:
                break

def umi_fastq_iter(path):
    """ Iterate through a fastaq file that is the result of MIGEC-assemble
        the UMI is in the fastq header
    """
    UMI_regex = re.compile('^@MIG UMI:(.*):([0-9]+)$')

    for header, seq, qual in fastq_iter(path):
        match = UMI_regex.match(header)
        umi, count = match.groups()
        yield umi, int(count), seq, qual

def clean_file(in_path, out_path, contaminated_umi):
    """ Given a seq of contaminated barcodes, filter all reads with one.
    """
    with open(out_path, 'w+') as out:
        # read_total, read_taken = 0, 0
        # clono_total, clono_taken = 0, 0
        for umi, count, seq, qual in umi_fastq_iter(in_path):
            # clono_total += 1
            # read_total += count
            if umi not in contaminated_umi:
                # clono_taken += 1
                # read_taken += count
                header = '@MIG UMI:%s:%i' %(umi, count)
                read = '%s\n%s\n+\n%s\n' % (header, seq, qual)
                out.write(read)

    # return clono_taken, clono_total, read_taken, read_total

def parse_config(path):
    result = []
    with open(path, 'r') as f:
        for line in f:
            split = line.strip().split('\t')
            if len(split) != 2:
                raise ValueError('Config file must contain two tab-seperated values per line')

            result.append(split)
    return result

def edit_dist(s1, s2):
    return Levenshtein.distance(s1, s2) / float(max(len(s1), len(s2)))

def write_arr_tsv(path, names, arr):
    with open(path, 'w+') as fout:

        # Write header
        fout.write('\t'+'\t'.join(names)+'\n')

        for name, row in zip(names, arr):
            row_str = '\t'.join(map(str, row))
            fout.write(name+'\t'+row_str+'\n')

def main(out_dir, files, t):

    if not os.path.exists(out_dir):
        print('Making directory:', out_dir)
        os.mkdir(out_dir)

    # UMI -> [sequences] . Map UMI to all sequences that have it.
    all_umi = defaultdict(list)

    # file_index -> set of UMI
    contaminated_umi = defaultdict(set)
    umi_counts = defaultdict(Counter)

    # Store the pairwise contamination by unique clonotypes and total reads.
    n = len(files)
    contamination_seqs = np.zeros((n, n), dtype='uint64')
    contamintation_reads = np.zeros((n, n), dtype='uint64')

    total_seqs = np.zeros(n, dtype='uint64')
    total_reads = np.zeros(n, dtype='uint64')

    ############################################################################
    print('Reading files...')

    for i, (name, f) in enumerate(files):
        print('\t', name)
        for umi, count, seq, _ in umi_fastq_iter(f):
            all_umi[umi].append((i, seq))
            umi_counts[i][umi] += count

            total_seqs[i] += 1
            total_reads[i] += count

    print('Finished Reading Files')
    print('There are %i UMI\'s' % len(all_umi))

    ############################################################################

    print('Calculating contaminated UMI\'s..')

    n_contam_pairings = 0

    for i, (umi, sequences) in enumerate(all_umi.items()):
        if i % 100000 == 0:
            print("{0:.0f}%".format(100 * i / float(len(all_umi))))

        if len(sequences) == 1:
            continue

        for (fi, s1), (fj, s2) in combinations(sequences, 2):

            if edit_dist(s1, s2) <= t:
                contaminated_umi[fi].add(umi)
                contaminated_umi[fj].add(umi)

                # Log information.
                contamination_seqs[fi, fj] += 1
                contamination_seqs[fj, fi] += 1

                contamintation_reads[fi, fj] += umi_counts[fi][umi]
                contamintation_reads[fj, fi] += umi_counts[fj][umi]

                n_contam_pairings += 1

    print('\nThere are %i contaminated UMI pairings!\n' % n_contam_pairings)

    ############################################################################
    print('Cleaning files...')

    for i, (name, file) in enumerate(files):
        out_path = os.path.join(out_dir, name+'.cleaned.fastq')
        clean_file(file, out_path, contaminated_umi[i])
        print('Wrote:', out_path)

    ############################################################################
    print('Writing result tables..')

    names = [name for name, _ in files]
    write_arr_tsv(os.path.join(out_dir, 'contamination_seqs.tsv'), names, \
                                                             contamination_seqs)
    write_arr_tsv(os.path.join(out_dir, 'contamination_reads.tsv'), names, \
                                                        contamintation_reads)

    write_arr_tsv(os.path.join(out_dir, 'contamination_seqs_fraction.tsv'), \
                    names, contamination_seqs.astype('float64')/total_seqs)
    write_arr_tsv(os.path.join(out_dir, 'contamination_reads_fraction.tsv'), \
                    names, contamintation_reads.astype('float64')/total_reads)

    print('Finished all.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', "--out", help="Directory for output files.")
    parser.add_argument('-c', "--config", help="Path to config file containing all fastq files.")
    parser.add_argument('-t', default=.05, help="Threshold for sequence similarity.")

    args = parser.parse_args()
    main(args.out, parse_config(args.config), args.t)
