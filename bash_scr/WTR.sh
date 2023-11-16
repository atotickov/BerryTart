#!/bin/bash

export PATH=/mnt/tank/scratch/skliver/common/mustelidae/atotik/variant_calling/samtools-1.18/bin:/mnt/tank/scratch/atotikov/tools/windowmasker-1.0.0:/nfs/home/atotikov/tools/Biocrutch/scripts/RepeatMasking:/nfs/home/atotikov/tools/BerryTart/python_scr:/mnt/tank/scratch/skliver/common/mustelidae/atotik/variant_calling/bedtools-2.31.0/bin:${PATH}

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
    
        echo $(date)" | Getting started on the ${sample_prefix} sample:"; echo;
    
        echo $(date)" | Step (1/5) | samtools faidx"
        samtools faidx ${sample}
        echo $(date)" | Step (1/5) | samtools faidx | Done"; echo;
    
        mkdir WindowMasker/
        cd WindowMasker/
    
        echo $(date)" | Step (2/5) | WindowMasker"
        windowmasker -in ../${sample} -infmt fasta -mk_counts -parse_seqids -out ${sample_prefix}.counts
        windowmasker -in ../${sample} -infmt fasta -ustat ${sample_prefix}.counts -outfmt interval -parse_seqids -out ${sample_prefix}.interval
        WindowMasker.py -i ${sample_prefix}.interval -o ${sample_prefix}.interval
        gzip ${sample_prefix}.counts
        gzip ${sample_prefix}.interval
        echo $(date)" | Step (2/5) | WindowMasker | Done"; echo;
    
        cd ../
        mkdir TRF/
        cd TRF/
        mv ../${sample} .
    
        echo $(date)" | Step (3/5) | TRF"
        
        # TRF застрянет в длинной центромерной области. Если нужно идентифицировать тандемный повтор, особенно в сборке T2T, нужно установить параметр -l с высоким значением, например -l 6
        # Параметры:
                # Match (matching weight): 2 
                # Mismatch (mismatching penalty): 7 
                # Delta (indel penalty): 7 
                # PM (match probability): 80 
                # PI (indel probability): 10 
                # Minscore (minimum alignment score): 50 
                # MaxPeriod (maximum period size): 2000 
                # -l (maximum TR length expected (in millions, eg, -l 3 or -l=3 for 3 million)): 10

        singularity run --bind $(pwd) --pwd $(pwd) ../../dfam-tetools-latest.sif trf ${sample} 2 7 7 80 10 50 2000 -l 10 -d -h
        mv ${sample}.2.7.7.80.10.50.2000.dat ${sample_prefix}.2.7.7.80.10.50.2000.l10.dat
        TRF.py -i ${sample_prefix}.2.7.7.80.10.50.2000.l10.dat -o ${sample_prefix}.2.7.7.80.10.50.2000.l10.dat
        gzip ${sample_prefix}.2.7.7.80.10.50.2000.l10.dat
        GFF_to_TAB.py ${sample_prefix}.2.7.7.80.10.50.2000.l10.dat.gff.gz ${sample_prefix}.2.7.7.80.10.50.2000.l10.dat.gff.tab
        echo $(date)" | Step (3/5) | TRF | Done"; echo;
    
        mv ${sample} ../
        cd ../
        mkdir RepeatMasker/
        cd RepeatMasker/
        mv ../${sample} .
    
        echo $(date)" | Step (4/5) | RepeatMasker"
        singularity run --bind $(pwd) --pwd $(pwd) ../../dfam-tetools-latest.sif RepeatMasker -pa 2 -species carnivora ${sample}
        mv ${sample}.out ${sample_prefix}.out
        RepeatMasker.py -i ${sample_prefix}.out -o ${sample_prefix}.out
        gzip ${sample}.masked
        gzip ${sample_prefix}.out
        echo $(date)" | Step (4/5) | RepeatMasker | Done"; echo;
    
        mv ${sample} ../
        cd ../
    
        echo $(date)" | Step (5/5) | Merging and sorting GFFs, masking repeats in ${sample} file"
        zcat TRF/${sample_prefix}.2.7.7.80.10.50.2000.l10.dat.gff.gz RepeatMasker/${sample_prefix}.out.gff.gz WindowMasker/${sample_prefix}.interval.gff.gz | sort -k1,1 -k4,4n -k5,5n | gzip > ${sample_prefix}.trf.repeatmasker.windowmasker.gff.gz
        bedtools maskfasta -soft -bed ${sample_prefix}.trf.repeatmasker.windowmasker.gff.gz -fi ${sample} -fo ${sample_prefix}.masked.fasta && gzip ${sample_prefix}.masked.fasta
        echo $(date)" | Step (5/5) | Merging and sorting GFFs, masking repeats in ${sample} file | Done"
    
        gzip ${sample}
    
        cd ../
done




# ------------------------------------------------------------------------
# https://github.com/Dfam-consortium/TETools
# conda create -n TETools -c conda-forge singularity
# conda activate TETools
# singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
# srun --cpus-per-task=32 --mem=256G --time=150:00:00 -w nd-1 --pty bash
# conda activate TETools
# bash WTR.sh -i 1.fasta 2.fasta 3.fasta 4.fasta 2>&1 | tee WTR.log
# ------------------------------------------------------------------------
