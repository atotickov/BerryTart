#!/bin/bash

export PATH=/home/skliver/atotikov/tools/gatk-4.4.0.0:${PATH}
export TOOLS="/nfs/home/atotikov/tools"

print_usage() {
        echo "Alternate Reference Maker (GATK):"
        echo "  -w      full path to genome assembly file (without '/')."
        echo "  -f      genome assembly file in .fasta format."
        echo "  -v      full path to vcf files."
}

assembly_path=""
assembly=""
vcf_path=""
vcfs=()

while getopts 'w:f:v:' flag; do
        case "${flag}" in
                w) assembly_path="${OPTARG}" ;;
                f) assembly="${OPTARG}" ;;
                b) vcf_path="${OPTARG}" ;;
                p) IFS=' ' read -ra vcfs <<< "${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done

workflow_path=$(pwd)





samtools faidx GCF_009829155.1_mMusErm1.Pri_genomic.fna

picard CreateSequenceDictionary R=GCF_009829155.1_mMusErm1.Pri_genomic.fna

/home/skliver/atotikov/tools/gatk-4.4.0.0/gatk IndexFeatureFile -I E26sub12.muserm.sub12.correct.filtered.masked.snp.vcf


/home/skliver/atotikov/tools/gatk-4.4.0.0/gatk FastaAlternateReferenceMaker --output E26sub12.AlternateReference.fasta --reference ../GCF_009829155.1_mMusErm1.Pri_genomic.fna --variant E26sub12.muserm.sub12.correct.filtered.masked.snp.vcf --use-iupac-sample E26sub12 --showHidden true



for i in "${vcfs[@]}"; do
        







mkdir ${workflow_path}/no_ChrX
mkdir ${workflow_path}/all_Chr && cd ${workflow_path}/all_Chr




