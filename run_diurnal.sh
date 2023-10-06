#!/bin/bash

path=`pwd`
cd $path'/scripts/'
pathdata=$path'/data/'

# case sens hour
echo $1
echo $2
echo $3

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

typ=$(IFS=, ; echo "${typs[*]}")
trac=$(IFS=, ; echo "${tracer[*]}")
python plot_stats.py $typ $trac $sens $hours $case $name $subcloud $cloud $sum
