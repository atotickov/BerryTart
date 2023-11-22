#!/usr/bin/env python3

# BLAST_with_colors.py -i blastn.GRCh38.consensuses.ev0.001.tab -o blastn.GRCh38.consensuses.ev0.001.color.tab -l legend.tab

import argparse

DEFAULT_COLORS = ['lightcoral', 'aquamarine', 'purple', 'sandybrown', 'gold', 'indianred', 'deepskyblue', 'red', 'lawngreen',
                  'deeppink', 'olive', 'chocolate', 'crimson', 'cadetblue', 'brown', 'lime', 'darkorange', 'dodgerblue',
                  'violet', 'darkgoldenrod', 'slategray', 'mediumslateblue', 'yellowgreen', 'orangered', 'darkturquoise',
                  'darkorchid', 'orange', 'cornflowerblue', 'palevioletred', 'coral']

def add_color_to_data(input_table, output_table, legend_table, colors_list):
    with open(input_table) as data_file:
        data = data_file.readlines()

    unique_clusters = {line.strip().split('\t')[0] for line in data}
    cluster_color_map = {cluster: color for cluster, color in zip(sorted(unique_clusters), colors_list)}

    with open(output_table, "w") as output_file:
        output_file.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tcolor\n")
        for line in data:
            cluster = line.strip().split('\t')[0]
            color = cluster_color_map[cluster]
            output_file.write(f"{line.strip()}\t{color}\n")

    with open(legend_table, "w") as legend_file:
        for cluster, color in cluster_color_map.items():
            legend_file.write(f"{color}\t{cluster}\n")

    return cluster_color_map

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add colors to data and create a legend')
    parser.add_argument('-i', '--input', dest='input_table', required=True, help='Input blastn 6 table')
    parser.add_argument('-o', '--output', dest='output_table', required=True, help='Output table name')
    parser.add_argument('-l', '--legend', dest='legend_table', required=True, help='Legend table name')
    parser.add_argument('--colors', nargs='*', help='List of colors')

    args = parser.parse_args()
    colors = args.colors if args.colors else DEFAULT_COLORS
    cluster_color_mapping = add_color_to_data(args.input_table, args.output_table, args.legend_table, colors)
