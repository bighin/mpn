#!/bin/bash

for i in *.ini ;
do
    OUTFILE=$(echo $i | sed "s/\.ini/\.log/g")
    sbatch --export=INIFILE=$i --output=$OUTFILE mpn.sbatch
    echo $i
done
