#!/bin/bash

export DYLD_LIBRARY_PATH=$(root-config --libdir):$DYLD_LIBRARY_PATH

for i in {0..35}; do
  echo "Running with argument $i"
  root -l -q -b "angular_dist_17O_half_dets.C($i)"
done
