#!/bin/bash

# Modifies output of using fragments with ropebwt3 into a 
# file equivalent to running ropebwt3 first without KeBaB.
# Usage: ./ropefix.sh input_file > output_file

if [ $# -ne 1 ]; then
    echo "Usage: $0 [MEM_FILE]"
    echo "[MEM_FILE] is the output of running ropebwt3 with KeBaB fragments"
    exit 1
fi

awk -F'\t' '{
    split($1, a, /[:-]/)
    o = a[2] - 1
    print a[1] "\t" $2 + o "\t" $3 + o "\t" $4
}' "$1"
