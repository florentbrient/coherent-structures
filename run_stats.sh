#!/bin/bash
#SBATCH --export=NONE
#SBATCH -N 1            # nodes number (=NBP)
#SBATCH -n 1            # CPUs number (on all nodes) (=NBP*TPN)
#SBATCH -t 04:00:00     # time limit

# Echo des commandes
#ulimit -c 0
#ulimit -s unlimited
# Arrete du job des la premiere erreur
#set -e
# Nom de la machine
#hostname

path=`pwd` 
cd $path'/scripts/'
pathdata=$path'/data/'

# case sens hour
echo $1
echo $2
echo $3
echo $4

case=$1
sens=$2
hours=$3
name='V0301'
svt=($4 $5 $6)

subcloud=0
cloud=0
sum='or'

typs=(updraft downdraft downdraft downdraft)
tracer=(${svt[0]}_WT ${svt[0]}_WT ${svt[1]}_WT ${svt[2]}_WT)
len=${#typs[@]}

hours=${hours//,/' '}
echo $hours

for i in $(seq 0 $len)
#for i in $(seq $len $len)
do
typ="${typs[$i]}"
trac="${tracer[$i]}"
if [ $i -eq $len ];
then
typ=$(IFS=, ; echo "${typs[*]}")
trac=$(IFS=, ; echo "${tracer[*]}")
fi
echo $typ
echo $trac

for hour in $hours #"${hours[@]}"
do
echo 'run'
python stats_flux.py $typ $trac $sens $hour $case $name $subcloud $cloud $sum
done
done


