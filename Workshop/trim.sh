#!/bin/bash

begin=`date +%s`
echo $HOSTNAME
printf ">polyA\nAAAAAAAAAAAAA\n>polyT\nTTTTTTTTTTTTT\n" | gzip - >  bbmap/resources/polyA.fa.gz
for sample in `cat samples.txt`
do
R1=${sample}_R1.fastq.gz

bbmap/bbduk.sh in=00-RawData/${R1} out=01-Trimmed/${R1} ref=bbmap/resources/polyA.fa.gz,bbmap/resources/truseq.fa.gz k=13 ktrim=r forcetrimleft=11 useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20

done

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

