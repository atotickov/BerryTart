#!/bin/bash

print_usage() {
        echo "  -l      Links to samples download from NCBI ("SRRlink1 SRRlink2 SRRlink3")."
}

SRR_links=();

while getopts 'l:' flag; do
        case "${flag}" in
                l) IFS=' ' read -ra SRR_links <<< "${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done

for link in "${SRR_links[@]}"; do
        SRR_ID=$(basename $link)
        echo $(date)" | SRR: ${SRR_ID} | Link: ${link} | Downloading";
        axel -n 100 ${link};
        echo $(date)" | SRR: ${SRR_ID} | Link: ${link} | Downloading | Done";

        echo $(date)" | SRR: ${SRR_ID} | Rename ${SRR_ID} -> ${SRR_ID}.sra";
        mv ${SRR_ID} ${SRR_ID}.sra;

        echo $(date)" | SRR: ${SRR_ID} | Splitting ${SRR_ID}.sra into fastq files";
        fastq-dump --split-3 ${SRR_ID}.sra 2> ${SRR_ID}.sra.log;
        #rm ${SRR_ID}.sra;
        echo $(date)" | SRR: ${SRR_ID} | Splitting ${SRR_ID}.sra into fastq files | Done";

        echo $(date)" | SRR: ${SRR_ID} | Gzipping fastq files";
        pigz -p 32 ${SRR_ID}_1.fastq;
        pigz -p 32 ${SRR_ID}_2.fastq;
        echo $(date)" | SRR: ${SRR_ID} | Gzipping fastq files | Done";
done




# ------------------------------------------------------------------------------------------------------------------------------------------------------------
#sbatch -w orthrus-1 SRA_slurm.sh

# Содержимое SRA_slurm.sh:

##!/bin/bash -i
##SBATCH --job-name=srrs                                                 # Job name
##SBATCH --mail-type=END                                                 # Mail events
##SBATCH --mail-user=a.totickov1@gmail.com                               # Where to send mail
##SBATCH --cpus-per-task=32                                              # Number of CPU cores per task (max 32)
##SBATCH --mem=256gb                                                     # Job memory request (max 256gb)
##SBATCH --time=150:00:00                                                # Time limit hrs:min:sec
##SBATCH --output=/logs/SRA_slurm.log
##SBATCH --error=/logs/SRA_slurm.err
#squeue; hostname; date;

#conda activate SRAtoolkit; # mamba create -n SRAtoolkit -c bioconda sra-tools; conda install -c conda-forge axel; conda install -c conda-forge pigz

#cd /path/; pwd;

#bash SRA.sh -l "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17072712/SRR17072712 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17072713/SRR17072713 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17072714/SRR17072714 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17072715/SRR17072715" 2>&1 | tee logs/SRA.info

#date;

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
