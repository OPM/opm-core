#!/bin/bash

source `dirname $0`/build-opm-core.sh

# Upstream revisions
declare -a upstreams
upstreams=(opm-parser
           opm-material
           opm-output)

declare -A upstreamRev
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-output]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

if grep -q "ert=" <<< $ghprbCommentBody
then
  ERT_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ert=([0-9]+).*/\1/g'`/merge
fi

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  OPM_COMMON_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

for upstream in $upstreams
do
  if grep -q "$upstream=" <<< $ghprbCommentBody
  then
    upstreamRev[$upstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*$upstream=([0-9]+).*/\1/g"`/merge
  fi
done

echo "Building with ert=$ERT_REVISION opm-common=$OPM_COMMON_REVISION opm-parser=${upstreamRev[opm-parser]} opm-material=${upstreamRev[opm-material]} opm-output=${upstreamRev[opm-output]} opm-core=$sha1"

build_opm_core
test $? -eq 0 || exit 1

# If no downstream builds we are done
if ! grep -q "with downstreams" <<< $ghprbCommentBody
then
  cp serial/build-opm-core/testoutput.xml .
  exit 0
fi

source $WORKSPACE/deps/opm-common/jenkins/setup-opm-data.sh
source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

# Downstream revisions
declare -a downstreams
downstreams=(opm-grid
             opm-simulators
             opm-upscaling
             ewoms)

declare -A downstreamRev
downstreamRev[opm-grid]=master
downstreamRev[opm-simulators]=master
downstreamRev[opm-upscaling]=master
downstreamRev[ewoms]=master

build_downstreams opm-core

test $? -eq 0 || exit 1
