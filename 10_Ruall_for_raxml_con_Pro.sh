#!/bin/bash

cd /mnt/lustre/users/wyunjin/leiming/Plant_phylogenomes/results/2015_1_27_raxml_trees/

echo "TIME START:"
date 

filename=$(echo $1 |awk -F "_" '{print $1"_"$2"_"$3"_"$4"_"$5}')
echo $filename
outfile_path=$(pwd)
echo  $outfile_path

mkdir $outfile_path/$filename\_pro
mkdir $outfile_path/$filename\_pro_con_modelselection

cd $outfile_path/$filename\_pro_con_modelselection 
ln -s ../$1    .

perl   /mnt/lustre/users/wyunjin/leiming/Plant_phylogenomes/results/2014_10_30_homologous/ProtienModelselection.pl  $1  >   $filename\_con_pro_model.txt 

con_model=$(awk '{if($0~/Best/){print $NF}}'  $filename\_con_pro_model.txt)
cd ../

~/leiming/src/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3   -f  a  -T 15  -x 12345  -p 12345 -m  PROTCAT$con_model  -#  200  -s  $1  -w  $outfile_path/$filename\_pro   -n  $filename\_pro.tr

echo "TIME END";
