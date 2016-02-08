#!/usr/bin/env bash
set -e

pushd . > /dev/null
git clone https://github.com/starseeker/SuperLU.git
cd SuperLU
mkdir build
cd build
cmake -D CMAKE_INSTALL_PREFIX=.. -D SUPERLU_BUILD_EXAMPLES=OFF -D SUPERLU_ENABLE_TESTING=OFF ../
make install
popd > /dev/null
