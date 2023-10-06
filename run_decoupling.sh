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

case='FIRE'
sens='Ls2x0'
svt=(SVT004 SVT005 SVT006)
hours=('003','006','009','012','015','018','021','024')
name='V0301'

hours=${hours//,/' '}
echo $hours


for hour in $hours #"${hours[@]}"
do
echo 'run'
python decoupling_index.py $sens $hour $case $name $subcloud
done


