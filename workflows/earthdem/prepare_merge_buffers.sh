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

echo "Linking '_reg_fill_merge.mat' files"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 00-matfiles/ --src-suffix _reg_fill.mat --dst-suffix _reg_fill_merge.mat --verbose
sleep_for_duration


echo "Creating neighbor index"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic create-neighbor-index "$UTM_ZONE"
sleep_for_duration

echo "Preforming dryrun for merge buffers by row"
conda run -n "$CONDA_ENV_NAME" --live-stream earthdem-mosaic merge-buffers "$UTM_ZONE" ./all_supertiles.txt row --slurm --dryrun

echo "Processing complete."
echo ""
echo "Run the following command to submit to the cluster:"
echo "    conda run -n $CONDA_ENV_NAME --live-stream earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt row --slurm"
exit 0