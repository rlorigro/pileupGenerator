#!/bin/bash
if [ $# -ne 8 ];
then
echo "usage: "$(basename $0) "[bam-file]" "[ref-file]" "[vcf-file]" "[confident-vcf]" "[coverage]" "[window-size]" "[output_dir]" "[max_threads]"
exit
fi

bam_file=$1
ref_file=$2
vcf_file=$3
confident_vcf_file=$4
coverage=$5
window_size=$6
output_dir=$7
max_threads=$8
now=$(date '+%d%m%Y%H%M%S')
output_dir=$output_dir'run-'$now
mkdir $output_dir
echo $output_dir
mkdir tmp
chrs=(1)
# chrs=(3)
for i in ${chrs[@]};
    do
    current_output_dir=$output_dir/'chr'$i/
    mkdir $current_output_dir
    echo "Starting chr" $i
    python3 pileupGenerator.py --bam $bam_file --ref $ref_file --vcf $vcf_file  --vcf_confident $confident_vcf_file --vcf_region $i --coverage $coverage --window_size $window_size --output_dir $current_output_dir --max_threads $max_threads 2>tmp/progress-$i.txt &
    done
