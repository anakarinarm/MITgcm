#!/bin/bash
set -e

cd ../build_lr
../../../tools/genmake2 -mods ../code_lr \
  -optfile ../../../tools/build_options/darwin_arm64_gfortran.opt -mpi
make depend
make -j 4
