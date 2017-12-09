#!/bin/bash
if [ $# -ne 7 ];
then
echo "usage: "$(basename $0) "[bam-file]" "[ref-file]" "[vcf-file]" "[coverage]" "[window-size]" "[output_dir]" "[max_threads]"
exit
fi

bam_file=$1
ref_file=$2
vcf_file=$3
coverage=$4
window_size=$5
output_dir=$6
max_threads=$7
now=$(date '+%d%m%Y%H%M%S')
output_dir=$output_dir'run-'$now
mkdir $output_dir
echo $output_dir
mkdir tmp
chrs=(3)
# chrs=(3)
for i in ${chrs[@]};
    do
    current_output_dir=$output_dir/'chr'$i/
    mkdir $current_output_dir
    echo "Starting chr" $i
    python3 pileupGenerator.py --bam $bam_file --ref $ref_file --vcf $vcf_file --vcf_region $i --coverage $coverage --window_size $window_size --output_dir $current_output_dir --max_threads $max_threads 2>tmp/progress-$i.txt &
    done
