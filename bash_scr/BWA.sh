# ------------------------------------------------------------------------------------------------------------------------------------------------------------
#sbatch -w orthrus-1 BWA_slurm.sh

# Содержимое BWA_slurm.sh:

##!/bin/bash -i
##SBATCH --job-name=musnig                                               # Job name
##SBATCH --mail-type=END                                                 # Mail events
##SBATCH --mail-user=a.totickov1@gmail.com                               # Where to send mail
##SBATCH --cpus-per-task=32                                              # Number of CPU cores per task (max 32)
##SBATCH --mem=256gb                                                     # Job memory request (max 256gb)
##SBATCH --time=150:00:00                                                # Time limit hrs:min:sec
##SBATCH --output=/mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musnig/alignment/new/logs/BWA_slurm.log
##SBATCH --error=/mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musnig/alignment/new/logs/BWA_slurm.err
#squeue; hostname; date;

#conda activate BWA; # bwa v.0.7.17-r1188, mosdepth v0.3.3

#cd /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musnig/alignment/new/; pwd;

#bash BWA.sh -a /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musnig/assembly/GCF_022355385.1_MUSNIG.SB6536_genomic.fna -s "SB7462sub12 SB6536sub12 SB8055sub12 SRR1508215sub12 SRR1508214sub12 SRR1508749sub12 SRR1508750sub12" -f /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musnig/alignment/new 2>&1 | tee logs/BWA.info

#date;

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

# Содержимое BWA.sh:


#!/bin/bash

export PATH=/mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/bcftools-1.18/bin:/mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/samtools-1.18/bin:${PATH};
export TOOLS="/nfs/home/atotikov/tools/";


print_usage() {
        echo "  -a      path to genome assembly file in .fasta format (with .len, .syn. .whitelist, .renamelist, and bwa index output)."
        echo "  -s      samples ID (example: "SB8055sub12 SB8156sub12 SB0109sub12")."
        echo "  -f      full path to paired sample files in .fastq.gz format (example: /path/to (without '/') SB8055sub12.filtered_1.fastq.gz, SB8055sub12.filtered_2.fastq.gz, etc.)."
}

assembly=""
samples=()
fq_path=""

while getopts 'a:s:f:l:w:' flag; do
        case "${flag}" in
                a) assembly="${OPTARG}" ;;
                s) IFS=' ' read -ra samples <<< "${OPTARG}" ;;
                f) fq_path="${OPTARG}" ;;
                l) assembly_len="${OPTARG}" ;;
                w) assembly_whitelist="${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done




for i in "${samples[@]}"; do
  echo $(date)" | Sample ${i} | (1/7) | Align paired-end reads to genome assembly, sorting and marking duplicates";
  bwa mem -t 25 ${assembly} ${fq_path}/${i}.filtered_1.fastq.gz ${fq_path}/${i}.filtered_2.fastq.gz -R "@RG\tID:${i}\tPU:x\tSM:${i}\tPL:Illumina\tLB:x" | samtools fixmate -@ 10 -m - - | samtools sort -@ 15 -m 12G | samtools markdup -@ 6 - ${i}.bam;
  
  echo $(date)" | Sample ${i} | (2/7) | Indexing of the created alignment file";
  samtools index ${i}.bam;
  
  echo $(date)" | Sample ${i} | (3/7) | Genome coverage calculation";
  mosdepth --threads 32 mosdepth.${i} ${i}.bam;
  
  echo $(date)" | Sample ${i} | (4/7) | Calculating whole-genome coverage statistics";
  python3 $TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i mosdepth.${i}.per-base.bed.gz --whole-genome-stats -o mosdepth.${i}
  
  echo $(date)" | Sample ${i} | (5/7) | Calculating nonoverlapping coverage statistics";
  python3 $TOOLS/Biocrutch/scripts/Coverage/coverage_statistics.py -i mosdepth.${i}.per-base.bed.gz --nonoverlapping-windows-stats -f 1000000 -o mosdepth.${i}
  
  echo $(date)" | Sample ${i} | (6/7) | Create mask file";
  python3 $TOOLS/MAVR/scripts/alignment/coverage/generate_mask_from_coverage_bed.py -c mosdepth.${i}.per-base.bed.gz -m $(cat mosdepth.${i}_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') -x 2.5 -n 0.33 -o mosdepth.${i}.max250.min33.bed

  echo $(date)" | Sample ${i} | (7/7) | Create mask file";
  python3 $TOOLS/MACE/scripts/draw_coverage.py --scaffold_column_name '#scaffold' --hide_track_label --window_column_name frame --coverage_column_name median -i mosdepth.${i}_1000000_windows_stats.csv -o ${i}.1Mb.track --subplots_adjust_left 0.12 --figure_width 12 -l "Sample ${i}: coverage (-w = 1mb)" -m $(cat mosdepth.${i}_whole_genome_stats.csv | sed -n 2p | awk '{print $2}') --window_size 1000000 --scaffold_length_file ${assembly}.len --scaffold_white_list ${assembly}.whitelist --scaffold_ordered_list ${assembly}.orderlist --scaffold_syn_file ${assembly}.syn --colormap jet --rounded;

done

