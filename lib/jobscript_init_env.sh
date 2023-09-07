#!/bin/bash

while (( $# != 0 )); do
    env_request="$1"
    env_load_cmd=''

    if [ "$env_request" == 'python' ]; then
        env_load_cmd="source ~/.bashrc ; conda activate pgc"

    elif [ "$env_request" == 'gdal' ]; then
        env_load_cmd="source ~/.bashrc ; conda activate pgc"

    elif [ "$env_request" == 'matlab' ]; then
#        env_load_cmd="module load matlab/2019a"
        env_load_cmd=""

    else
        echo "Unsupported environment request: '${env_request}'"
        exit 1
    fi

    if [ -n "$env_load_cmd" ]; then
        echo "Loading '${env_request}' environment with the following command(s):"
        printf '%s\n' "$env_load_cmd"
        eval "$env_load_cmd"
        rc=$?
        if (( rc != 0 )); then
            echo "WARNING: Non-zero return code (${rc}) from env load command(s)"
        fi
        echo
    fi

    shift
done
