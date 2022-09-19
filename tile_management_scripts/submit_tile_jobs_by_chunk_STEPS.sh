#!/bin/bash

big_tilelist="tilelists/tilelists_trex_altered/earthdem_tiles_intersecting_trex_priorities_mar2022.txt"
chunk_size=500
next_chunk_jobnum_thresh=100
jobname_prefix='t2t_'

tilelist_chunk_base="${big_tilelist/.txt/_}"

# Split large tilelist into smaller chunks
split --numeric-suffixes --lines "$chunk_size" "$big_tilelist" "$tilelist_chunk_base"

# That produced files like:
# ${tilelist_chunk_base}_00
# ${tilelist_chunk_base}_01
# ...

# Rename them to be nice
for f in ${tilelist_chunk_base}_* ; do mv "$f" "${f}.txt" ; done

# Submit tiles by chunk, moving to the next chunk when there are fewer than ${next_chunk_jobnum_thresh} jobs in the queue
for tilelist in ${tilelist_chunk_base}_*.txt; do
    python ~/common_repos/setsm_postprocessing_pgc/batch_tiles2tif_v4.py results/output_tiles/ "$tilelist" 2 earthdem --pbs --hold
    while true; do
        sleep 5m
        njobs=$(qstat -u "$USER" | grep -c "$jobname_prefix")
        if (( njobs < next_chunk_jobnum_thresh )); then
            break
        fi
        echo "Number of '${jobname_prefix}' jobs = ${njobs}, will submit next chunk when < ${next_chunk_jobnum_thresh}";
    done
done
