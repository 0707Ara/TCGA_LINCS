#!/bin/bash

ncol_origin=""
nrow_origin=""
ncol_sub=""
nrow_sub=""
arr=""
ls
cd TCGA
for file in *; do
	cd $file
		for f in *; do
			if [[ $f == *.gz ]]; then
				ncol_origin=`zcat $f | head -n 1 | tr $'\t' '\n' | wc -l`
                                nrow_origin=`zcat $f | wc -l`
                                arr+=("$f")
                                arr+=("$ncol_origin")
                                arr+=("$nrow_origin")
			elif [[ $f == *.txt ]]; then
				ncol_sub=`awk '{print NF}' $f | sort -nu | tail -n 1`
				nrow_sub=`wc -l < $f`
				arr+=("$f")
				arr+=("$ncol_sub")
				arr+=("$nrow_sub")
			fi
		done

	cd ../
done

for value in "${arr[@]}";do
	echo $value >>/home/arajo0707/TCGA_LINCS/temp.txt
done
