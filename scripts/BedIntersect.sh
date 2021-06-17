#!/usr/bin/env bash

bedtools sort -i $1 > $2
bedtools sort -i $3 > $4
bedtools intersect -a $2 -b $4 -loj -wa -sorted >$5 #no strandness required!
