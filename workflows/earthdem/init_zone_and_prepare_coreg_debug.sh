#!/usr/bin/env bash

set -euo pipefail

readonly CONDA_ENV_NAME="earthdem-mosaic"
readonly SLEEP_DURATION="10s"

error() {
  local line=$1      # Line number from ${LINENO}
  local cmd=$2       # Command from $BASH_COMMAND
  echo "Error on line $line: Command '$cmd' failed with status $?" >&2
}
trap 'error ${LINENO} "$BASH_COMMAND"' ERR

sleep_for_duration() {
    echo "Sleeping for $SLEEP_DURATION seconds before proceeding..."
    sleep "$SLEEP_DURATION"
}

if [[ -v UTM_ZONE ]]; then
    echo "Processing zone: $UTM_ZONE"
else
    echo "Environment variable UTM_ZONE not set"
    exit 1
fi

if [[ -v EARTHDEM_MOSAIC_ENV_FILE ]]; then
    echo "Using config: $EARTHDEM_MOSAIC_ENV_FILE"
    conda run -n "$CONDA_ENV_NAME" earthdem-mosaic show-settings
else
    echo "Environment variable EARTHDEM_MOSAIC_ENV_FILE not set"
    exit 1
fi

sleep_for_duration

echo "Creating working directories"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic create-working-dirs "$UTM_ZONE" --verbose
sleep_for_duration


echo "Linking source matfiles"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic link-source-matfiles "$UTM_ZONE" --verbose
sleep_for_duration

echo "Creating supertile list"
find "./$UTM_ZONE/00-matfiles/" -maxdepth 1 -type d -name "utm*" | sed "s|./$UTM_ZONE/00-matfiles/||" | sort > "./$UTM_ZONE/all_supertiles.txt"
sleep_for_duration

echo "Linking matfiles to 10-coregistration-debug directory"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic link-files-to-stage --src "./$UTM_ZONE/00-matfiles" --dst "./$UTM_ZONE/10-coregistration-debug" --src-suffix ".mat" --verbose

echo "Linking fin files to 10-coregistration-debug directory"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic link-files-to-stage --src "./$UTM_ZONE/00-matfiles" --dst "./$UTM_ZONE/10-coregistration-debug" --src-suffix ".fin" --verbose
sleep_for_duration

echo "Preforming dryrun for coreg-debug stage"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic coreg-debug "$UTM_ZONE" "./$UTM_ZONE/all_supertiles.txt" --slurm --dryrun

echo "Processing complete. Run the following command to submit to the cluster:"
echo "cd $UTM_ZONE && conda run -n $CONDA_ENV_NAME --live-stream earthdem-mosaic coreg-debug $UTM_ZONE ./all_supertiles.txt --slurm"
exit 0