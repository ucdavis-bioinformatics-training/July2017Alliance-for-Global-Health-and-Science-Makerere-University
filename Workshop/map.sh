#!/bin/bash

begin=`date +%s`
echo $HOSTNAME
mkdir 02-Salmon
# index the transcriptome
Salmon-0.8.2_linux_x86_64/bin/salmon index --type quasi -t References/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa -i References/Saccharomyces_cerevisiae.R64-1-1.index

for sample in `cat samples.txt`
do
R1=${sample}_R1.fastq.gz

Salmon-0.8.2_linux_x86_64/bin/salmon quant -i References/Saccharomyces_cerevisiae.R64-1-1.index --noLengthCorrection -l SF -r 01-Trimmed/${R1} -o 02-Salmon/${R1}.quant

done

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

