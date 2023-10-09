#!/bin/bash

# ------------------------------------------------------------------------
# https://github.com/Dfam-consortium/TETools
# conda create -n TETools -c conda-forge singularity pigz
# conda activate TETools
# singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
# srun --cpus-per-task=32 --mem=256G --time=150:00:00 -w nd-1 --pty bash
# conda activate TETools
# bash WTR.sh -i 1.fasta 2.fasta 3.fasta 4.fasta 2>&1 | tee WTR.log
# ------------------------------------------------------------------------

export PATH=/path/to/samtools-1.18/bin:/path/to/windowmasker-1.0.0:/path/to/Biocrutch/scripts/RepeatMasking:/path/to/bedtools-2.31.0/bin:${PATH}

if [ "$1" != "-i" ]; then
    echo "Usage: $0 -i 1.fasta 2.fasta 3.fasta ..."
    exit 1
fi

shift 1

for sample in "$@"; do
    if [ ! -f "$sample" ]; then
        echo "'$sample' file not found, moving on to the next file."
        continue
    fi
    
    sample_prefix="${sample%.*}"
    mkdir "$sample_prefix"
    cd "$sample_prefix"
    
    mv ../${sample} .
    
    echo $(date)" | Getting started on the ${sample_prefix} sample:"
    
    echo ' '
    echo $(date)" | Step (1/5) | samtools faidx"
    samtools faidx ${sample}
    echo $(date)" | Step (1/5) | samtools faidx | Done"
    
    echo ' '
    mkdir WindowMasker/
    cd WindowMasker/
    
    echo $(date)" | Step (2/5) | WindowMasker"
    windowmasker -in ../${sample} -infmt fasta -mk_counts -parse_seqids -out ${sample_prefix}.counts
    windowmasker -in ../${sample} -infmt fasta -ustat ${sample_prefix}.counts -outfmt interval -parse_seqids -out ${sample_prefix}.interval
    WindowMasker.py -i ${sample_prefix}.interval -o ${sample_prefix}.interval
    pigz -p 32 ${sample_prefix}.counts
    pigz -p 32 ${sample_prefix}.interval
    echo $(date)" | Step (2/5) | WindowMasker | Done"
    
    echo ' '
    cd ../
    mkdir TRF/
    cd TRF/
    mv ../${sample} .
    
    echo $(date)" | Step (3/5) | TRF"
    singularity run --bind $(pwd) --pwd $(pwd) ../../dfam-tetools-latest.sif trf ${sample} 2 7 7 80 10 50 500 -d -h
    mv ${sample}.2.7.7.80.10.50.500.dat ${sample_prefix}.2.7.7.80.10.50.500.dat
    TRF.py -i ${sample_prefix}.2.7.7.80.10.50.500.dat -o ${sample_prefix}.2.7.7.80.10.50.500.dat
    pigz -p 32 ${sample_prefix}.2.7.7.80.10.50.500.dat
    echo $(date)" | Step (3/5) | TRF | Done"
    
    echo ' '
    mv ${sample} ../
    cd ../
    mkdir RepeatMasker/
    cd RepeatMasker/
    mv ../${sample} .
    
    echo $(date)" | Step (4/5) | RepeatMasker"
    singularity run --bind $(pwd) --pwd $(pwd) ../../dfam-tetools-latest.sif RepeatMasker -pa 2 -species carnivora ${sample}
    mv ${sample}.out ${sample_prefix}.out
    RepeatMasker.py -i ${sample_prefix}.out -o ${sample_prefix}.out
    pigz -p 32 ${sample}.masked
    pigz -p 32 ${sample_prefix}.out
    echo $(date)" | Step (4/5) | RepeatMasker | Done"
    
    echo ' '
    mv ${sample} ../
    cd ../
    
    echo $(date)" | Step (5/5) | Merging and sorting GFFs, masking repeats in ${sample} file"
    zcat TRF/${sample_prefix}.2.7.7.80.10.50.500.dat.gff.gz RepeatMasker/${sample_prefix}.out.gff.gz WindowMasker/${sample_prefix}.interval.gff.gz | sort -k1,1 -k4,4n -k5,5n | gzip > ${sample_prefix}.trf.repeatmasker.windowmasker.gff.gz
    bedtools maskfasta -soft -bed ${sample_prefix}.trf.repeatmasker.windowmasker.gff.gz -fi ${sample} -fo ${sample_prefix}.masked.fasta && gzip ${sample_prefix}.masked.fasta
    echo $(date)" | Step (5/5) | Merging and sorting GFFs, masking repeats in ${sample} file | Done"
    
    pigz -p 32 ${sample}
    
    cd ../

done

