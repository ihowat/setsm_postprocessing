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


echo "Building 20-no-slope-filter VRT"
cd 20-no-slope-filter/
gdalbuildvrt ./20-no-slope-filter_browse.vrt $(find -type f -name "*_browse.tif" | paste -sd " ")
cd ..

echo "Building 30-yes-slope-filter VRT"
cd 30-yes-slope-filter/
gdalbuildvrt ./30-yes-slope-filter_browse.vrt $(find -type f -name "*_browse.tif" | paste -sd " ")
cd ..

echo "Building slope filter review GeoPackage"
mamba run -n "$MAMBA_ENV_NAME" --live-stream earthdem-mosaic slope-filter-review "$UTM_ZONE" --verbose

echo "Processing complete."