#!/bin/bash

export PATH=/mnt/tank/scratch/atotikov/tools/psmc:/mnt/tank/scratch/atotikov/tools/psmc/utils:${PATH}
export TOOLS="/nfs/home/atotikov/tools"

print_usage() {
        echo "Correct variant calling:"
        echo "  -f      genome assembly file in .fasta format (with .fasta file there should be a .fasta.fai file)."
        echo "  -b      alignment file in .bam format (with .bam file there should be a .bam.bai file)."
        echo "  -s      sample ID."
        echo "  -m      sample whole genome stats (mosdepth.SRR6963888.sub12_whole_genome_stats.csv)."
}

while getopts 'f:b:s:m:' flag; do
        case "${flag}" in
                f) assembly="${OPTARG}" ;;
                b) alignment="${OPTARG}" ;;
                s) prefix="${OPTARG}" ;;
                m) median="${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done

echo  $(date)" | Stage 1 | Variant Calling";

mkdir -p split/ split//bcf/;

python3 $TOOLS/MAVR/scripts/sequence/prepare_region_list.py -r ${assembly}.fai -s -m 1500000 -n 1 -g samtools -x 1000 2>/dev/null | parallel -j 2 "samtools mpileup -C 50 -uf ${assembly} -r {} ${alignment} | bcftools view -b -c - > split//bcf//tmp.{#}.bcf";

echo  $(date)" | Stage 2 | Bcftools cat - concatenation";

bcftools cat `ls split//bcf//tmp.*.bcf | sort -V` >> ${prefix}.bcf;

rm -r split/;

echo  $(date)" | Stage 3 | Bcftools view - bcf to vcf conversation";

bcftools view ${prefix}.bcf | pigz -p 32 >> ${prefix}.vcf;

echo  $(date)" | Stage 4 | Consensus file preparation";

zcat ${prefix}.vcf.gz | vcfutils.pl vcf2fq -d $(cat ${median} | sed -n 2p | awk '{print $2 / 3}') -D $(cat ${median} | sed -n 2p | awk '{print $2 * 2.5}') | pigz -p 32 > ${prefix}.fq.gz;

echo  $(date)" | Stage 5 | Fasta-like consensus file preparation";

fq2psmcfa -q20 ${prefix}.fq.gz > ${prefix}.diploid.psmcfa;

echo  $(date)" | Stage 6 | Fasta-like consensus file preparation for bootstrapping";

splitfa ${prefix}.diploid.psmcfa > ${prefix}.diploid.split.psmcfa;

echo  $(date)" | Stage 7 | General";

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${prefix}.diploid.psmc ${prefix}.diploid.psmcfa;

echo  $(date)" | Stage 8 | Bootstrapping";

seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc ${prefix}.diploid.split.psmcfa > seq100.txt;
parallel -j 32 < seq100.txt;
cat round-*.psmc > ${prefix}.round.psmc;
rm round-*;
rm seq100.txt;

echo  $(date)" | Stage 9 | Generation time calculation"

mkdir G2 && cd G2/;

psmc_plot.pl -u2.2e-9 -g2 -R -p diploid ../${prefix}.diploid.psmc;

mv diploid.0.txt ${prefix}.u2_2e-9.g2.diploid.txt;
rm diploid.*;

psmc_plot.pl -u2.2e-9 -g2 -R -p round ../${prefix}.round.psmc;

for i in round.*.txt; do cat $i >> ${prefix}.u2_2e-9.g2.round.txt; done

rm round.*;
cd ../;


#conda activate portugal;

#bash test.sh -f /mnt/tank/scratch/skliver/common/mustelidae/mustela_erminea/genome/assemblies/GCF_009829155.1_mMusErm1.Pri_genomic.fna -b /mnt/tank/scratch/skliver/common/mustelidae/mustela_erminea/genome/alignment/12/SRR6963888.sub12.bam -s SRR6963888 -m /mnt/tank/scratch/skliver/common/mustelidae/mustela_erminea/genome/alignment/12/mosdepth.SRR6963888.sub12_whole_genome_stats.csv
