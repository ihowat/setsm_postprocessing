#!/usr/bin/env bash

set -euo pipefail

readonly SLEEP_DURATION="5s"

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
    earthdem-mosaic show-settings
else
    echo "Environment variable EARTHDEM_MOSAIC_ENV_FILE not set"
    exit 1
fi


echo "Copying '_reg_fill_merge.mat' files"
sleep_for_duration
earthdem-mosaic link-files-to-stage --src 00-matfiles/ --dst 00-matfiles/ \
                                    --src-suffix _reg_fill.mat --dst-suffix _reg_fill_merge.mat \
                                    --copy --verbose


echo "Creating neighbor index"
sleep_for_duration
earthdem-mosaic create-neighbor-index "$UTM_ZONE"


echo "Preforming dryrun for merge buffers by row"
sleep_for_duration
earthdem-mosaic merge-buffers "$UTM_ZONE" ./all_supertiles.txt row --slurm --dryrun


echo "Processing complete."
echo ""
echo "Run the following command to submit to the cluster:"
echo "    earthdem-mosaic merge-buffers $UTM_ZONE ./all_supertiles.txt row --slurm"
exit 0