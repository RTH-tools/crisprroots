#!/bin/bash

input=$1;
VAR=$(find $input/* -type d | wc -l)
find "$@"  -type f -name "*.bam" -printf '%f\n' | sed 's/^[^.]*.//g'  | sort | uniq -c | awk -v VAR="$VAR" '{if ($1 == VAR) print $2}'
