#!/bin/sh
(bwa mem -t 12  -R '@RG\tID:foo\tSM:A2_B2\tLB:library1' $1 $2 $3 | samtools sort  --write-index -T $4 -o $5 -)  2> $6