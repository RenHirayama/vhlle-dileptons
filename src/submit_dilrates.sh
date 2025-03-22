#!/bin/bash

path=$1
energy=$2
evmin=$3
evmax=$4

for evnum in $(seq $evmin $evmax);do
  sed  -e "s|HYDROPATH|${path}|g" -e "s/ENERGY/${energy}/g"  -e "s/EVNUM/${evnum}/g" \
       template_dilrates.sh > job_dilrate_en"${energy}"-${evnum}.sh
done
