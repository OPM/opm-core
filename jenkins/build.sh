#!/bin/bash

source `dirname $0`/build-opm-core.sh

# Upstream revisions
declare -a upstreams
upstreams=(opm-parser
           opm-material)

declare -A upstreamRev
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

build_opm_core
test $? -eq 0 || exit 1

cp serial/build-opm-core/testoutput.xml .
