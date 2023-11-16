bash test.sh -f /mnt/tank/scratch/skliver/common/mustelidae/atotik/PSMC/muserm/assembly/GCF_009829155.1_mMusErm1.Pri_genomic.fna -B "SRR6963880sub12" -b /mnt/tank/scratch/skliver/common/mustelidae/atotik/PSMC/muserm/alignment -s "muserm" -M "mosdepth.SRR6963880sub12" -m /mnt/tank/scratch/skliver/common/mustelidae/atotik/PSMC/muserm/alignment -g 1 -u 2.2e-9 -c "NC_045635.1" 2>&1 | tee test.info


bash PSMC.sh -f /mnt/tank/scratch/skliver/common/mustelidae/atotik/PSMC/muserm/assembly/GCF_009829155.1_mMusErm1.Pri_genomic.fna -B "SRR6963880sub12 SRR6963883sub12 SRR6963884sub12 SRR6963885sub12 SRR6963886sub12 SRR6963887sub12 SRR6963888sub12 SRR6963889sub12 E26sub12 T136sub12" -b /mnt/tank/scratch/skliver/common/mustelidae/atotik/PSMC/muserm/alignment -s "muserm" -M "mosdepth.SRR6963880sub12 mosdepth.SRR6963883sub12 mosdepth.SRR6963884sub12 mosdepth.SRR6963885sub12 mosdepth.SRR6963886sub12 mosdepth.SRR6963887sub12 mosdepth.SRR6963888sub12 mosdepth.SRR6963889sub12 mosdepth.E26sub12 mosdepth.T136sub12" -m /mnt/tank/scratch/skliver/common/mustelidae/atotik/PSMC/muserm/alignment -g 1 -u 2.2e-9 -c "NC_045635.1" 2>&1 | tee PSMC.info

#!/bin/bash

export PATH="/mnt/tank/scratch/atotikov/tools/psmc:/mnt/tank/scratch/atotikov/tools/psmc/utils:${PATH}"
export TOOLS="/nfs/home/atotikov/tools"

print_usage() {
    echo "Correct variant calling:"
    echo "  -f      Genome assembly file in .fasta format (with .fasta file there should be a .fasta.fai file)."
    echo "  -B      Alignment id names (with .bam files there should be a .bam.bai files; example: 'SB8055sub12 SB8156sub12 SB0109sub12')."
    echo "  -b      Full path to alignment files in .bam format (without '/'; with .bam files there should be a .bam.bai files)."
    echo "  -s      Species Latin name."
    echo "  -M      Samples whole genome stat files id (mosdepth.SRR6963888.sub12)."
    echo "  -m      Full path to samples whole genome stat files (without '/'; mosdepth.SRR6963888.sub12_whole_genome_stats.csv)."
    echo "  -g      Generation time."
    echo "  -u      Mutation rate."
    echo "  -c      ChrX scaffold name."
}

assembly=""
bam=()
bams_path=""
prefix=""
median=()
medians_path=""
generation_time=""
mutation_rate=""
ChrX_ID=""

while getopts 'f:B:b:s:M:m:g:u:c:' flag; do
        case "${flag}" in
                f) assembly="${OPTARG}" ;;
                B) IFS=' ' read -ra bam <<< "${OPTARG}" ;;
                b) bams_path="${OPTARG}" ;;
                s) prefix="${OPTARG}" ;;
                M) IFS=' ' read -ra median <<< "${OPTARG}" ;;
                m) medians_path="${OPTARG}" ;;
                g) generation_time="${OPTARG}" ;;
                u) mutation_rate="${OPTARG}" ;;
                c) ChrX_ID="${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done

workflow_path=$(pwd)
mkdir ${workflow_path}/no_ChrX
mkdir ${workflow_path}/all_Chr && cd ${workflow_path}/all_Chr

for ((i=0; i<${#bam[@]}; i++)); do

    sample=${bam[i]%sub12}
    mkdir ${workflow_path}/all_Chr/${sample} && cd ${workflow_path}/all_Chr/${sample}

    echo  $(date)" | Sample: ${sample}"
    echo  $(date)" | Stage 1 | Variant Calling";

    mkdir -p ${workflow_path}/all_Chr/${sample}/split/ ${workflow_path}/all_Chr/${sample}/split//bcf/

    python3 $TOOLS/MAVR/scripts/sequence/prepare_region_list.py -r ${assembly}.fai -s -m 1500000 -n 1 -g samtools -x 1000 2>/dev/null | parallel -j 32 "samtools mpileup -C 50 -uf ${assembly} -r {} ${bams_path}/${bam[i]}.bam | bcftools view -b -c - > split//bcf//tmp.{#}.bcf";

    echo  $(date)" | Stage 1 | Variant Calling | Done";echo;
    echo  $(date)" | Stage 2 | Concatenation";

    bcftools cat `ls split//bcf//tmp.*.bcf | sort -V` >> ${prefix}.${sample}.bcf;

    rm -r ${workflow_path}/all_Chr/${sample}/split/;
    echo  $(date)" | Stage 2 | Concatenation | Done";echo;
    echo  $(date)" | Stage 3 | BCF to VCF conversation";

    bcftools view ${prefix}.${sample}.bcf | gzip > ${prefix}.${sample}.vcf.gz;

    rm ${prefix}.${sample}.bcf
    echo  $(date)" | Stage 3 | BCF to VCF conversation | Done";
    echo  $(date)" | Stage 4 | Consensus file preparation | Stat file: ${median[i]}_whole_genome_stats.csv | -d: $(awk -v value=$(cat ${medians_path}/${median[i]}_whole_genome_stats.csv | sed -n 2p | awk '{print $2 / 3}') 'BEGIN { printf "%.0f\n", value }') | -D: $(awk -v value=$(cat ${medians_path}/${median[i]}_whole_genome_stats.csv | sed -n 2p | awk '{print $2 * 2.5}') 'BEGIN { printf "%.0f\n", value }')";

    zcat ${prefix}.${sample}.vcf.gz | vcfutils.pl vcf2fq -d $(awk -v value=$(cat ${medians_path}/${median[i]}_whole_genome_stats.csv | sed -n 2p | awk '{print $2 / 3}') 'BEGIN { printf "%.0f\n", value }') -D $(awk -v value=$(cat ${medians_path}/${median[i]}_whole_genome_stats.csv | sed -n 2p | awk '{print $2 * 2.5}') 'BEGIN { printf "%.0f\n", value }') | gzip > ${prefix}.${sample}.fq.gz;

    echo  $(date)" | Stage 4 | Consensus file preparation | Done";
    echo  $(date)" | Stage 5 | Fasta-like consensus file preparation";

    fq2psmcfa -q20 ${prefix}.${sample}.fq.gz > ${prefix}.${sample}.diploid.psmcfa;

    echo  $(date)" | Stage 5 | Fasta-like consensus file preparation | Done";
    echo  $(date)" | Stage 6 | Fasta-like consensus file preparation for bootstrapping";

    splitfa ${prefix}.${sample}.diploid.psmcfa > ${prefix}.${sample}.diploid.split.psmcfa;

    echo  $(date)" | Stage 6 | Fasta-like consensus file preparation for bootstrapping | Done";
    echo  $(date)" | Stage 7 | General population history calculating";

    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${prefix}.${sample}.diploid.psmc ${prefix}.${sample}.diploid.psmcfa;

    echo  $(date)" | Stage 7 | General population history calculating | Done";
    echo  $(date)" | Stage 8 | 100 bootstrapping";

    seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc ${prefix}.${sample}.diploid.split.psmcfa > seq100.txt;
    parallel -j 32 < seq100.txt;
    cat round-*.psmc > ${prefix}.${sample}.round.psmc;
    rm round-*.psmc;
    rm seq100.txt;

    echo  $(date)" | Stage 8 | 100 bootstrapping | Done";
    echo  $(date)" | Stage 9 | Population history calculation | Generation time: ${generation_time} | Mutation rate: ${mutation_rate}"

    mkdir ${workflow_path}/all_Chr/${sample}/G${generation_time}_U${mutation_rate} && cd ${workflow_path}/all_Chr/${sample}/G${generation_time}_U${mutation_rate};

    psmc_plot.pl -u ${mutation_rate} -g ${generation_time} -R diploid ${workflow_path}/all_Chr/${sample}/${prefix}.${sample}.diploid.psmc;

    mv diploid.0.txt ${prefix}.${sample}.G${generation_time}_U${mutation_rate}.diploid.txt;
    rm diploid.*;

    psmc_plot.pl -u ${mutation_rate} -g ${generation_time} -R round ${workflow_path}/all_Chr/${sample}/${prefix}.${sample}.round.psmc;

    for r in round.*.txt; do cat $r >> ${prefix}.${sample}.G${generation_time}_U${mutation_rate}.round.txt; done
    rm round.*;
    echo  $(date)" | Stage 9 | Population history calculation | Generation time: ${generation_time} | Mutation rate: ${mutation_rate} | Done"

    mkdir ${workflow_path}/no_ChrX/${sample} && cd ${workflow_path}/no_ChrX/${sample}

    echo  $(date)" | Stage 10 | Excluding ${ChrX_ID} scaffold"

    zcat ${workflow_path}/all_Chr/${sample}/${prefix}.${sample}.fq.gz | paste - - - - | grep -v -F -e "$ChrX_ID" | tr "\t" "\n" | gzip > ${prefix}.${sample}.no_ChrX.fq.gz

    echo  $(date)" | Stage 10 | Excluding ${ChrX_ID} scaffold | Done"
    echo  $(date)" | Stage 11 | Fasta-like consensus file preparation";

    fq2psmcfa -q20 ${prefix}.${sample}.no_ChrX.fq.gz > ${prefix}.${sample}.no_ChrX.diploid.psmcfa;

    echo  $(date)" | Stage 11 | Fasta-like consensus file preparation | Done";
    echo  $(date)" | Stage 12 | Fasta-like consensus file preparation for bootstrapping";

    splitfa ${prefix}.${sample}.no_ChrX.diploid.psmcfa > ${prefix}.${sample}.no_ChrX.diploid.split.psmcfa;

    echo  $(date)" | Stage 12 | Fasta-like consensus file preparation for bootstrapping | Done";
    echo  $(date)" | Stage 13 | General population history calculating";

    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${prefix}.${sample}.no_ChrX.diploid.psmc ${prefix}.${sample}.no_ChrX.diploid.psmcfa;

    echo  $(date)" | Stage 13 | General population history calculating | Done";
    echo  $(date)" | Stage 14 | 100 bootstrapping";

    seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc ${prefix}.${sample}.no_ChrX.diploid.split.psmcfa > seq100.txt;
    parallel -j 32 < seq100.txt;
    cat round-*.psmc > ${prefix}.${sample}.no_ChrX.round.psmc;
    rm round-*.psmc;
    rm seq100.txt;

    echo  $(date)" | Stage 14 | 100 bootstrapping | Done";
    echo  $(date)" | Stage 15 | Population history calculation | Generation time: ${generation_time} | Mutation rate: ${mutation_rate}"

    mkdir ${workflow_path}/no_ChrX/${sample}/G${generation_time}_U${mutation_rate} && cd ${workflow_path}/no_ChrX/${sample}/G${generation_time}_U${mutation_rate};

    psmc_plot.pl -u ${mutation_rate} -g ${generation_time} -R diploid ${workflow_path}/no_ChrX/${sample}/${prefix}.${sample}.no_ChrX.diploid.psmc;

    mv diploid.0.txt ${prefix}.${sample}.G${generation_time}_U${mutation_rate}.no_ChrX.diploid.txt;
    rm diploid.*;

    psmc_plot.pl -u ${mutation_rate} -g ${generation_time} -R round ${workflow_path}/no_ChrX/${sample}/${prefix}.${sample}.no_ChrX.round.psmc;

    for r in round.*.txt; do cat $r >> ${prefix}.${sample}.G${generation_time}_U${mutation_rate}.no_ChrX.round.txt; done
    rm round.*;
    echo  $(date)" | Stage 15 | Population history calculation | Generation time: ${generation_time} | Mutation rate: ${mutation_rate} | Done"

    cd ${workflow_path}/all_Chr/;
done























