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
				#extract number of rows and columns from .gz file and write to nrow.txt file
				ncol_origin=`zcat $f | head -n 1 | tr $'\t' '\n' | wc -l`
                                nrow_origin=`zcat $f | wc -l`
                                echo -n $f >>/home/arajo0707/TCGA_LINCS/nrow.txt
                                echo -e -n ' \t ' >>/home/arajo0707/TCGA_LINCS/nrow.txt
				echo -n $ncol_origin>>/home/arajo0707/TCGA_LINCS/nrow.txt
                                echo -e -n ' \t '>>/home/arajo0707/TCGA_LINCS/nrow.txt
				echo -n $nrow_origin>>/home/arajo0707/TCGA_LINCS/nrow.txt
				echo -e -n ' \t '>>/home/arajo0707/TCGA_LINCS/nrow.txt
			elif [[ $f == *.txt ]]; then
				#extract number of rows and columns from .txt file and write to nrow.txt file
				ncol_sub=`awk '{print NF}' $f | sort -nu | tail -n 1`
				nrow_sub=`wc -l < $f`
				echo -n $f >>/home/arajo0707/TCGA_LINCS/nrow.txt
                                echo -e -n ' \t ' >>/home/arajo0707/TCGA_LINCS/nrow.txt
                                echo -n $ncol_sub>>/home/arajo0707/TCGA_LINCS/nrow.txt
                                echo -e -n ' \t ' >>/home/arajo0707/TCGA_LINCS/nrow.txt
                                echo -n $nrow_sub>>/home/arajo0707/TCGA_LINCS/nrow.txt
				echo -e -n ' \t ' >>/home/arajo0707/TCGA_LINCS/nrow.txt
			fi
		done
	cd ../
done
