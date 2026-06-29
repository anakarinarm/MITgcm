#!/bin/bash
set -e

cd ../build_lr_obeta
../../../tools/genmake2 -mods ../code_lr_obeta \
  -optfile ../../../tools/build_options/darwin_arm64_gfortran.opt -mpi
make depend
make -j 4
