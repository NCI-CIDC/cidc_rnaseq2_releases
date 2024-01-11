#!/bin/bash

## Output the number of sequences 
size=$(cat $1 | awk 'FNR == 3 {{print}}' | grep -o '[[:digit:]]*')
size=$(($size/1000000))
## If size is greater than 10 million, downsample 20% of sequences
if [ $size -lt 10 ]
then
    downsampling_size=$size'M'
else
    downsampling_size=$(($size*2/10))'M'
fi
echo $downsampling_size > $2


VALUE=$(echo $downsampling_size | sed -r 's/M//g')
FACTOR=$(samtools idxstats $3 | cut -f3 | awk -v COUNT=$((VALUE*1000000)) 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]
then
   echo '[ERROR]: Requested number of reads exceeds total read count in' $3 '-- exiting' && exit 1
fi

## Subsample the bam by the calculated factor and sort the output bam
samtools view -@ $6 -s $FACTOR -b $3 | samtools sort -@ $6 -o $4

## Create index for the downsampled bam
samtools index -@ $6 $4 > $5
