#!/usr/bin/env python3

# python3 BLAST_to_HITS.py -i blastn.GRCh38.consensuses.ev0.001.tab -o blastn.GRCh38.consensuses.ev0.001.hits.tab

import argparse
from collections import defaultdict

def process_table(input_file, output_file):
    hits = defaultdict(int)
    scaffolds = {}
    clusters = {}
    all_scaffolds = []
    all_clusters = []

    with open(input_file, 'r') as file:
        for line in file:
            data = line.strip().split('\t')
            cluster = data[0]
            scaffold = data[1]
            hits[(cluster, scaffold)] += 1

            if scaffold not in scaffolds:
                scaffolds[scaffold] = 1
                all_scaffolds.append(scaffold)

            if cluster not in clusters:
                clusters[cluster] = 1
                all_clusters.append(cluster)

    all_scaffolds = sorted(all_scaffolds)

    with open(output_file, 'w') as output:
        output.write('\t' + '\t'.join(all_scaffolds) + '\n')
        for cluster in all_clusters:
            output.write(cluster + '\t')
            for scaffold in all_scaffolds:
                output.write(str(hits[(cluster, scaffold)]) + '\t')
            output.write('\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process hits table')
    parser.add_argument('-i', '--input', help='Input blastn 6 table file', required=True)
    parser.add_argument('-o', '--output', help='Output file for the new hits table', required=True)
    args = parser.parse_args()

    process_table(args.input, args.output)
