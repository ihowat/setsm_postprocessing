#!/bin/bash

if (( $# != 1 )); then
    echo "Script expects 1 argument: path to text file tilelist (each tile name listed on a separate line)"
    exit 1
fi

tilelist="$1"

base10() { printf '%s' "$((10#$1))"; }

while IFS= read -r tile || [ -n "$tile" ]; do
    row_col=$(echo "$tile" | grep -Eo '[0-9]{2}_[0-9]{2}')
    utm_prefix="${tile/${row_col}/}"
    row=$(base10 "$(echo "$row_col" | cut -d"_" -f1)")
    col=$(base10 "$(echo "$row_col" | cut -d"_" -f2)")
    for i in $((row-1)) $row $((row+1)); do
        for j in $((col-1)) $col $((col+1)); do
            if { (( i == row )) || (( j == col )); } && { (( i > 0 )) && (( j > 0 )); }; then
                tile_adj=$(printf '%s%02d_%02d' "$utm_prefix" $i $j )
                echo "$tile_adj"
            fi
        done
    done
done < "$tilelist" | sort -u
