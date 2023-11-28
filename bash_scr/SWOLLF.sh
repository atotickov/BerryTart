#!/bin/bash

print_usage() {
        echo "  -i      Assemblies file in fasta format ("assembly1.fasta assembly2.fasta assembly3.fasta")."
        echo "  -n      1n for each assembly ("22 21 16")" 
}

assemblies=();
ChrN=();

while getopts 'i:n:' flag; do
        case "${flag}" in
                i) IFS=' ' read -ra assemblies <<< "${OPTARG}" ;;
                n) IFS=' ' read -ra ChrN <<< "${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done

workflow_path=$(pwd)
echo "Workflow Path: ${workflow_path}"

for ((i=0; i<${#assemblies[@]}; i++)); do

        prefix=${assemblies[i]%.*}
        echo "Assembly: ${assemblies[i]} | Chr number (1n): ${ChrN[i]} | Prefix: ${prefix}"
        
        mkdir -p ${workflow_path}/${prefix}
        mv ${assemblies[i]} ${workflow_path}/${prefix}/
        cd ${workflow_path}/${prefix}/
        
        samtools faidx ${assemblies[i]}
        
        cat ${assemblies[i]}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 > ${prefix}.len 
        
        cat ${assemblies[i]}.fai | awk '{print $1"\t"$2}' | sort -nr -k2 | head -n ${ChrN[i]} > ${prefix}.lengths 
        
        cat ${prefix}.lengths | awk '{print $1}' > ${prefix}.whitelist 
        
        cat ${prefix}.whitelist | awk '{print $1"\t"$1}' > ${prefix}.syn 
        
        cat ${prefix}.syn | awk '{print $2}' > ${prefix}.orderlist
        
        cd ${workflow_path}/
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
