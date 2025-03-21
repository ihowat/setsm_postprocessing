#!/usr/bin/env bash

set -euo pipefail

readonly MAMBA_ENV_NAME="earthdem-mosaic"
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
    mamba run -n "$MAMBA_ENV_NAME" earthdem-mosaic show-settings
else
    echo "Environment variable EARTHDEM_MOSAIC_ENV_FILE not set"
    exit 1
fi

sleep_for_duration

echo "Linking .mat and .fin files to 20-no-slope-filter"
mamba run -n "$MAMBA_ENV_NAME" --live-stream earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 20-no-slope-filter/ --src-suffix _reg_fill_merge.mat --dst-suffix .mat
mamba run -n "$MAMBA_ENV_NAME" --live-stream earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 20-no-slope-filter/ --src-suffix .fin

echo "Preforming dryrun for no-slope-filter export"
mamba run -n "$MAMBA_ENV_NAME" --live-stream earthdem-mosaic export-final-tifs "$UTM_ZONE" ./all_supertiles.txt --slurm --dryrun
sleep_for_duration

echo "Linking .mat and .fin files to 30-yes-slope-filter"
mamba run -n "$MAMBA_ENV_NAME" --live-stream earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 30-yes-slope-filter/ --src-suffix _reg_fill_merge.mat --dst-suffix .mat
mamba run -n "$MAMBA_ENV_NAME" --live-stream earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 30-yes-slope-filter/ --src-suffix .fin

echo "Preforming dryrun for yes-slope-filter export"
mamba run -n "$MAMBA_ENV_NAME" --live-stream earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --apply-slope-filter --slurm --dryrun

echo "Processing complete."
echo ""
echo "Run the following commands to submit to the cluster:"
echo "    mamba run -n $MAMBA_ENV_NAME --live-stream earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --slurm"
echo "    mamba run -n $MAMBA_ENV_NAME --live-stream earthdem-mosaic export-final-tifs $UTM_ZONE ./all_supertiles.txt --slurm --apply-slope-filter"
exit 0