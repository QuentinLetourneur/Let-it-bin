#!/bin/bash

abun_file=$1

nb_col=$2

for i in `seq 3 $nb_col`;do
	file_name=`head -n 1 $abun_file | cut -f $i -d$'\t'`
	cut -f 1,$i $abun_file | tail -n +2 > $file_name
done
