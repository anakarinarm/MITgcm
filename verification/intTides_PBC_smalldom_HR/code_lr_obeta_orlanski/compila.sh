#!/bin/bash
set -e

cd ../build_lr_obeta_orlanski
../../../tools/genmake2 -mods ../code_lr_obeta_orlanski \
  -optfile ../../../tools/build_options/darwin_arm64_gfortran.opt -mpi
make depend
make -j 4
