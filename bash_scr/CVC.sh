#!/bin/bash

export PATH=/mnt/tank/scratch/skliver/common/mustelidae/atotik/variant_calling/bcftools-1.18/bin:/mnt/tank/scratch/skliver/common/mustelidae/atotik/variant_calling/samtools-1.18/bin:/mnt/tank/scratch/skliver/common/mustelidae/atotik/variant_calling/bedtools-2.31.0/bin:/mnt/tank/scratch/skliver/common/mustelidae/atotik/variant_calling/vcftools-0.1.16/bin:${PATH}
export TMPDIR=/mnt/tank/scratch/skliver/common/mustelidae/atotik/;
export TOOLS="/nfs/home/atotikov/tools"

print_usage() {
        echo "Correct variant calling:"
        echo "  -w      full path to genome assembly .syn .whitelist .orderedlist .len files (without '/')."
        echo "  -f      genome assembly file in .fasta format (with .fasta file there should be a .fasta.fai file)."
        echo "  -b      full path to alignment files in .bam format (with .bam file there should be a .bam.bai file)."
        echo "  -p      ploidy file."
        echo "  -s      sample file."
        echo "  -v      output vcf file prefix (example, musput.sub12.correct)."
        echo "  -m      full path to mask files in .bed format (without '/')."
        echo "  -o      species name ('Mustela putorius')."
}

while getopts 'w:f:b:p:s:v:m:o:' flag; do
        case "${flag}" in
                w) assembly_swol="${OPTARG}" ;;
                f) assembly="${OPTARG}" ;;
                b) alignment="${OPTARG}" ;;
                p) ploidy_file="${OPTARG}" ;;
                s) sample_file="${OPTARG}" ;;
                v) prefix_vcf="${OPTARG}" ;;
                m) mask="${OPTARG}" ;;
                o) species="${OPTARG}" ;;
                *) print_usage
                        exit 1 ;;
        esac
done




echo $(date)" | Stage 1 | Correct variant calling and concatination into one ${prefix_vcf}.vcf.gz file";

mkdir -p split/ split//mpileup/ split//bcf/;

echo "Step 1/2 | Correct variant calling";
python3 $TOOLS/MAVR/scripts/sequence/prepare_region_list.py -r ${assembly}.fai -s -m 1500000 -n 100 -g samtools -x 1000 2>/dev/null | parallel -j 32 "bcftools mpileup -d 250 -q 30 -Q 30 --adjust-MQ 50 -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP,SCR,INFO/SCR -O u -f ${assembly} -r {} `find ${alignment} -name *.bam -type f | tr "\n" " "` | bcftools call --ploidy-file ${ploidy_file} --samples-file ${sample_file} --group-samples - -m -O u -v -f GQ,GP > split//bcf//tmp.{#}.bcf";

echo "Step 2/2 | Concatination into one ${prefix_vcf}.vcf.gz file";
bcftools concat -O u --threads 32 `ls split//bcf//tmp.*.bcf | sort -V` | bcftools view -O z -o ${prefix_vcf}.vcf.gz -;

rm -r split/;
echo $(date)" | Stage 1 | Correct variant calling and concatination into one ${prefix_vcf}.vcf.gz file | Done"; echo;




echo $(date)" | Stage 2 | Common ${prefix_vcf}.vcf.gz file filtration"
mkdir ${prefix_vcf}.bcftools_filtration/ && cd ${prefix_vcf}.bcftools_filtration/
echo "! Filtration params: 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)'"

bcftools filter --threads 32 -S . -O z -o ${prefix_vcf}.filtered.vcf.gz --exclude 'QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)' ../${prefix_vcf}.vcf.gz

echo $(date)" | Stage 2 | Common ${prefix_vcf}.vcf.gz file filtration | Done"; echo;




echo $(date)" | Stage 3 | Separation of samples from the ${prefix_vcf}.filtered.vcf.gz file and masking"

touch samples.filtered.masked.txt

for sampleID in `bcftools query -l ${prefix_vcf}.filtered.vcf.gz`; do
        echo "Step 1/2 | Sample ${sampleID} separation";
        bcftools view --threads 32 --min-ac 1 --with-header -s ${sampleID} -O z -o ${sampleID}.${prefix_vcf}.filtered.vcf.gz ${prefix_vcf}.filtered.vcf.gz;

        echo "Step 2/2 | Sample ${sampleID} masking";
        bedtools intersect -header -v -a ${sampleID}.${prefix_vcf}.filtered.vcf.gz -b ${mask}/mosdepth.${sampleID}.max250.min33.bed > ${sampleID}.${prefix_vcf}.filtered.masked.vcf;
        echo ${sampleID}.${prefix_vcf}.filtered.masked.vcf >> samples.filtered.masked.txt;
done

echo $(date)" | Stage 3 | Separation of samples from the ${prefix_vcf}.filtered.vcf.gz file and masking | Done"; echo;




echo $(date)" | Stage 4 | Separation of genetic variants"

touch samples.filtered.masked.all.txt

for filtered_masked in $(cat samples.filtered.masked.txt); do
        echo "Sample ${filtered_masked%%.*}"

        echo "Step 1/6 | Indel separation";
        bcftools filter --threads 32 -i 'TYPE="indel"' -O z -o ${filtered_masked%.*}.indel.vcf.gz ${filtered_masked};
        echo ${filtered_masked%.*}.indel.vcf.gz >> samples.filtered.masked.all.txt;

        echo "Step 2/6 | Snp separation";
        bcftools filter --threads 32 -i 'TYPE="snp"' -O z -o ${filtered_masked%.*}.snp.vcf.gz ${filtered_masked};
        echo ${filtered_masked%.*}.snp.vcf.gz >> samples.filtered.masked.all.txt;

        echo "Step 3/6 | Hetero indel variants separation"
        bcftools filter --threads 32 -i 'FMT/GT="het"' -O z -o ${filtered_masked%.*}.indel.hetero.vcf.gz ${filtered_masked%.*}.indel.vcf.gz;
        echo ${filtered_masked%.*}.indel.hetero.vcf.gz >> samples.filtered.masked.all.txt;

        echo "Step 4/6 | Hetero snp variants separation"
        bcftools filter --threads 32 -i 'FMT/GT="het"' -O z -o ${filtered_masked%.*}.snp.hetero.vcf.gz ${filtered_masked%.*}.snp.vcf.gz;
        echo ${filtered_masked%.*}.snp.hetero.vcf.gz >> samples.filtered.masked.all.txt;

        echo "Step 5/6 | Homo indel variants separation"
        bcftools filter --threads 32 -i 'FMT/GT="hom"' -O z -o ${filtered_masked%.*}.indel.homo.vcf.gz ${filtered_masked%.*}.indel.vcf.gz;
        echo ${filtered_masked%.*}.indel.homo.vcf.gz >> samples.filtered.masked.all.txt;

        echo "Step 6/6 | Homo snp variants separation"
        bcftools filter --threads 32 -i 'FMT/GT="hom"' -O z -o ${filtered_masked%.*}.snp.homo.vcf.gz ${filtered_masked%.*}.snp.vcf.gz;
        echo ${filtered_masked%.*}.snp.homo.vcf.gz >> samples.filtered.masked.all.txt;

        echo "Sample ${filtered_masked%%.*} | Done"
done

rm samples.filtered.masked.txt
echo $(date)" | Stage 4 | Separation of genetic variants | Done"; echo;




echo  $(date)" | Stage 5 | Statistics calculation"

echo "Filtration type 1 ('QUAL < 20.0 || (FORMAT/SP > 60.0 | FORMAT/DP < 5.0 | FORMAT/GQ < 20.0)')" > stats.csv
echo "File Name,Number of variants,Unique genotypes" >> stats.csv

while read -r filename; do
        numvar=$(zcat "$filename" | grep -vP "^#" | cut -f 10 | sed 's/:.*//' | grep -v "\./\." | wc -l)
        unique=$(zcat "$filename" | grep -vP "^#" | cut -f 10 | sed 's/:.*//' | sort | uniq)
        echo ${filename},${numvar},${unique} >> stats.csv
done < samples.filtered.masked.all.txt

rm samples.filtered.masked.all.txt
echo  $(date)" | Stage 5 | Statistics calculation | Done"; echo;





echo  $(date)" | Stage 6 | Visualization of genetic variants in heat maps"

for file in *.snp.hetero.vcf.gz; do
        echo "File ${file} | Hetero SNPs | 100kb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.100kb.hetero --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Heterozygous SNPs (-w, -s = 100kb)" -w 100000 -s 100000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;

        echo "File ${file} | Hetero SNPs | 1mb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.1mb.hetero --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Heterozygous SNPs (-w, -s = 1mb)" -w 1000000 -s 1000000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;
done

for file in *.snp.homo.vcf.gz; do
        echo "File ${file} | Homo SNPs | 100kb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.100kb.homo --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Homozygous SNPs (-w, -s = 100kb)" -w 100000 -s 100000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;

        echo "File ${file} | Homo SNPs | 1mb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).snps.1mb.homo --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Homozygous SNPs (-w, -s = 1mb)" -w 1000000 -s 1000000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;
done

for file in *.indel.hetero.vcf.gz; do
        echo "File ${file} | Hetero Indels | 100kb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.100kb.hetero --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Heterozygous indels (-w, -s = 100kb)" -w 100000 -s 100000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;

        echo "File ${file} | Hetero Indels | 1mb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.1mb.hetero --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Heterozygous indels (-w, -s = 1mb)" -w 1000000 -s 1000000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;
done

for file in *.indel.homo.vcf.gz; do
        echo "File ${file} | Homo Indels | 100kb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.100kb.homo --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Homozygous indels (-w, -s = 100kb)" -w 100000 -s 100000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;

        echo "File ${file} | Homo Indels | 1mb";
        python3 $TOOLS/MACE/scripts/draw_variant_window_densities.py -i ${file} -o $(echo $file | cut -d'.' -f 1).indels.1mb.homo --density_thresholds 0.00,0.10,0.50,1.00,1.50,2.00,2.50,3.00,4.00,6.00 -l "Sample $(echo $file | cut -d'.' -f 1): Homozygous indels (-w, -s = 1mb)" -w 1000000 -s 1000000 -a ${assembly_swol}/*.whitelist -z ${assembly_swol}/*.orderlist -n ${assembly_swol}/*.len --scaffold_syn_file ${assembly_swol}/*.syn --syn_file_key_column 0 --syn_file_value_column 1 --colormap jet --hide_track_label --rounded --subplots_adjust_left 0.12;
done

echo  $(date)" | Stage 6 | Visualization of genetic variants in heat maps | Done"; echo;




echo  $(date)" | Stage 7 | Visualization of genetic variants density in violin plot"

tsv_folder=$(pwd)
output_file="${prefix_vcf}.violinplot.1mb.tab"
touch "$output_file"
tab=$'\t'
add_entry_to_violinplot() {
        local file_name="$1"
        local identifier="$2"
        local gender="$3"
        local full_path_with_file_name="${tsv_folder}/${file_name}"
        echo "$identifier ($gender)$tab$full_path_with_file_name" >> "$output_file"
}


new_files=("$tsv_folder"/*.snps.1mb.hetero.variant_counts.tsv)

for new_file in "${new_files[@]}"; do
        file_name=$(basename "$new_file")
        identifier=$(echo "$file_name" | cut -d'.' -f1)
        gender=$(awk -v id="$identifier" '$1 == id {print $2}' "$sample_file")
        if [ -n "$gender" ]; then
                add_entry_to_violinplot "$file_name" "$identifier" "$gender"
        fi
done

python3 $TOOLS/Biocrutch/scripts/Visualization/draw_violinplots.py -i ${prefix_vcf}.violinplot.1mb.tab -o ${prefix_vcf}.violinplot.1mb.tab -w 1000000 --figure_height 9 --yticklist 0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7 --title "${species}" --ylabel "Гетерозиготные SNP/тыс п.н." --figure_width_per_sample 0.7 --rotation 70 --ymin 0 --ymax 7 --font-size 14 --figure_grid

echo  $(date)" | Stage 7 | Visualization of genetic variants density in violin plot | Done"

# Визуализация рохов -- нужно сделать
# На файлах features.bed, которые получаются после отрисовки гетерозиготности:
# python3 $TOOLS/Biocrutch/scripts/ROH/get_ROH_regions.py -i SAMPLE.snp.hetero.features.bed -o SAMPLE.snp.hetero.features.roh
# далее визуализация на хромосомах через draw_features.py
# потом отрисовка кумулятивных графиков, ноут есть в папке юпитер ноутбук на компе (нужно настроить)




# ------------------------------------------------------------------------------------------------------------------------------------------------------------
#sbatch -w orthrus-1 CVC_slurm.sh

# Содержимое CVC_slurm.sh:

##!/bin/bash -i
##SBATCH --job-name=musput                                               # Job name
##SBATCH --mail-type=END                                                 # Mail events
##SBATCH --mail-user=a.totickov1@gmail.com                               # Where to send mail
##SBATCH --cpus-per-task=32                                              # Number of CPU cores per task (max 32)
##SBATCH --mem=256gb                                                     # Job memory request (max 256gb)
##SBATCH --time=150:00:00                                                # Time limit hrs:min:sec
##SBATCH --output=/path/to/varcall/logs/CVC_slurm.log
##SBATCH --error=/path/to/varcall/logs/CVC_slurm.err
#squeue; hostname; date;

#conda activate varcall

#cd /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/varcall/; pwd;

#bash CVC.sh -w /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/assembly -f /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/assembly/mustela_putorius.ragtag.fasta -b /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/alignment/ -p /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/varcall/ploidy.sub12.file -s /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/varcall/sample.sub12.file -v musput.sub12.correct -m /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/varcall/TESTS/masks -o 'Mustela putorius' 2>&1 | tee logs/CVC.info

#bash CVC_violinplot.sh -s /mnt/tank/scratch/skliver/common/mustelidae/atotik/vc/musput/varcall/sample.sub12.file -v musput.sub12.correct -o 'Mustela putorius' 2>&1 | tee logs/CVC_violinplot.info

#date;

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
