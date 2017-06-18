#!/bin/bash

cd /mnt/lustre/users/wyunjin/leiming/Plant_phylogenomes/results/2015_1_27_orthology_assignment_subtree
echo "TIME START"
date

filename=$(echo $1 |awk -F "_" '{print $1"_"$2"_"$3}')
echo $filename

outfile_path=$(pwd)
echo  $outfile_path

mkdir $outfile_path/$filename\_dna

~/leiming/src/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3  -f a  -T 12  -x 12345 -p 12345 -m  GTRCATI -#  100   -s  $1 -w  $outfile_path/$filename\_dna  -n  $filename\_gene.tr

echo "TIME END"
date

#__END__
