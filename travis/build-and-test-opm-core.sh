#!/usr/bin/env bash
set -e

pushd . > /dev/null
opm-core/travis/build-opm-core.sh
cd opm-core/build
ctest --output-on-failure
popd > /dev/null
