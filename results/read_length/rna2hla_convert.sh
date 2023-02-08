#!/bin/bash
# declare an array called array and define 3 vales
array=( 36 41 46 51 56 61 66 71 76 81 86 91 96 101 106 111 116 121 126 )
for i in "${array[@]}"
do
	cat rna2hla_read_length.csv | grep _"$i" >> phlat_"$i".csv 
	#echo "$i"
done
