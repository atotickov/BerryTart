#!/bin/bash

#bash PAR.sh -i "mosdepth.SB6536sub12.per-base.bed.gz mosdepth.SB8055sub12.per-base.bed.gz mosdepth.SRR1508214sub12.per-base.bed.gz mosdepth.SRR1508215sub12.per-base.bed.gz mosdepth.SRR1508749sub12.per-base.bed.gz" -s NC_081575.1 -w 10000 2>&1 | tee PAR.info
export TOOLS="/nfs/home/atotikov/tools/"

print_usage() {
        echo "  -i      mosdepth.*.mosdepth.per-base.bed.gz files (example: "mosdepth.SB6536sub12.per-base.bed.gz mosdepth.SB8055sub12.per-base.bed.gz mosdepth.SRR1508214sub12.per-base.bed.gz")."
        echo "  -s      chrX ID."
        echo "  -w      window size (usually: 10000)."
}

mosdepth_bedgz=()
chrX_ID=""
window_size=""

while getopts "i:s:w:" flag; do
        case "${flag}" in
                i) IFS=' ' read -ra mosdepth_bedgz <<< "${OPTARG}" ;;
                s) chrX_ID="${OPTARG}" ;;
                w) window_size="${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done

PWD=$(pwd)

for mbedgz in "${mosdepth_bedgz[@]}"; do
  echo $(date)" | File ${mbedgz} | (1/3) | Getting the sex chromosome data";

  # mosdepth.sampleID.per-base.bed.gz -> mosdepth.sampleID + _PAR
  mkdir ${PWD}/${mbedgz%.*.*.*}_PAR/;
  cd ${PWD}/${mbedgz%.*.*.*}_PAR/;

  zcat ../${mbedgz} | awk '{ if ($1 == "'${chrX_ID}'") print $0}' | gzip --stdout > ${mbedgz%.*.*}.chrX.bed.gz;
  
  echo $(date)" | File ${mbedgz} | (2/3) | Window stats";
  python3 $TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i ${mbedgz%.*.*}.chrX.bed.gz --tool-name mosdepth -n -f ${window_size} -o ${mbedgz%.*.*}.chrX;
  # remove header
  sed -i '1,1d' ${mbedgz%.*.*}.chrX_${window_size}_windows_stats.csv;
  
  echo $(date)" | File ${mbedgz} | (3/3) | Computation of pseudoautosomal region coordinates";
  python3 $TOOLS/Biocrutch/scripts/PAR/pseudoautosomal_region.py -f ${window_size} -i ${mbedgz%.*.*}.chrX_${window_size}_windows_stats.csv -s ${chrX_ID} -m $(cat ../${mbedgz%.*.*.*}_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') -o ${mbedgz%.*.*}.chrX | tee ${mbedgz%.*.*}.chrX_pseudo.log

  cd ../

done
