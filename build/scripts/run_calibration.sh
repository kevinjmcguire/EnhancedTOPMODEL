#!/bin/sh
# Pre-calibrate all gauges (2 at a time); the dryland basins get a longer search and two seeds.
# Run from build/:  sh scripts/run_calibration.sh      (~15-20 min on 2 cores)
cd "$(dirname "$0")/.."
{
  for b in 02017500 02028500 02038850 01669000 01013500 12040500 11143000 11266500 09352900; do echo "$b 1"; done
  for b in 08202700 09512280; do echo "$b 1"; echo "$b 2"; done
} | xargs -P 2 -n 2 sh -c 'node scripts/calibrate.js $0 $1 > /dev/null && echo "done $0 seed $1"'
