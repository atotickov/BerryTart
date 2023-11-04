#!/bin/bash

print_usage() {
    echo "  -i      sample id ('1 2 3 4 5')."
    echo "  -s      sampling value (10000000 by default)"
}

IDs=()
sampling="10000000"

while getopts 'i:s:' flag; do
    case "${flag}" in
        i) IFS=' ' read -ra IDs <<< "${OPTARG}" ;;
        s) sampling="${OPTARG}" ;;
        *) print_usage
           exit 1 ;;
    esac
done

for ID in "${IDs[@]}"; do
    echo "$(date) | Sampling ${ID}"

    seed=$(od -An -N2 -i /dev/urandom | awk '{print $1}')
    echo "$(date) | Random seed value for the sample paired reads with ${ID} ID: ${seed}"
    echo "$(date) | Sampling value: ${sampling}"

    input_1="${ID}.filtered_1.fastq.gz"
    input_2="${ID}.filtered_2.fastq.gz"

    output_1="${ID}.filtered_1.sub${sampling}.fastq.gz"
    output_2="${ID}.filtered_2.sub${sampling}.fastq.gz"

    seqtk sample -s${seed} ${input_1} ${sampling} | pigz -p 32 > ${output_1}
    seqtk sample -s${seed} ${input_2} ${sampling} | pigz -p 32 > ${output_2}

    echo "$(date) | The number of sequences (1) and bases (2) in the ${output_1} file:"
    seqtk size ${output_1}
    echo "$(date) | The number of sequences (1) and bases (2) in the ${output_2} file:"
    seqtk size ${output_2}
done



# conda activate SeqTK
# STK.sh -i "1 2 3 4" -s 10000000
