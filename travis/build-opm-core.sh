#!/usr/bin/env bash
set -e

pushd . > /dev/null
cd opm-core
mkdir build
cd build
cmake -D SUPERLU_ROOT=../../SuperLU ../
make
popd > /dev/null
