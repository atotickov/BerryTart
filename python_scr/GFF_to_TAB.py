#!/usr/bin/env python3

import argparse
import gzip

def read_gff(gff_file):
    print(f"Reading {gff_file}: ", end="")
    gff_data = []
    with gzip.open(gff_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 9:
                chr_id, source, feature, start, end, score, strand, frame, attributes = parts
                gff_data.append([chr_id, source, feature, start, end, score, strand, frame, attributes])
            else:
                print(f"Skipping invalid line: {line}")
    
    print(f"found {len(gff_data)} rows")
    return gff_data

def parse_attributes(attributes_str, attrsep=";", fieldsep="="):
    attr_dict = {}
    for attr in attributes_str.split(attrsep):
        key, value = attr.split(fieldsep, 1)
        key = key.strip()
        value = value.strip()
        attr_dict[key] = value
    return attr_dict

def gff_to_tab(gff_file, output_file, attrsep=";", fieldsep="="):
    gff_data = read_gff(gff_file)
    with open(output_file, 'w') as output:
        output.write("\t".join(["chr_id", "source", "feature", "start", "end", "score", "strand", "frame"]))

    with open(output_file, 'a') as output:
        output.write("\tID\tperiod\tcopies\tconsensus_size\tperc_match\tperc_indels\talign_score\tperc_A\tperc_C\tperc_G\tperc_T\tentropy\tcons_seq\trepeat_seq\n")
        for row in gff_data:
            attributes = parse_attributes(row[8], attrsep=attrsep, fieldsep=fieldsep)
            output.write("\t".join(row[:8]) + f"\t{attributes.get('ID', '')}\t{attributes.get('period', '')}\t{attributes.get('copies', '')}\t{attributes.get('consensus_size', '')}\t{attributes.get('perc_match', '')}\t{attributes.get('perc_indels', '')}\t{attributes.get('align_score', '')}\t{attributes.get('perc_A', '')}\t{attributes.get('perc_C', '')}\t{attributes.get('perc_G', '')}\t{attributes.get('perc_T', '')}\t{attributes.get('entropy', '')}\t{attributes.get('cons_seq', '')}\t{attributes.get('repeat_seq', '')}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert gzipped GFF file to tab-delimited format")
    parser.add_argument("input_file", help="Input gzipped GFF file")
    parser.add_argument("output_file", help="Output tab-delimited file")
    args = parser.parse_args()
    
    gff_to_tab(args.input_file, args.output_file)
