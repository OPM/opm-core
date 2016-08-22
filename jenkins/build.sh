#!/bin/bash

source `dirname $0`/build-opm-core.sh

# Upstream revisions
declare -a upstreams
upstreams=(ert
           opm-parser
           opm-material
           opm-output)

declare -A upstreamRev
upstreamRev[ert]=master
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-output]=master

OPM_COMMON_REVISION=master

build_opm_core
test $? -eq 0 || exit 1

cp serial/build-opm-core/testoutput.xml .
