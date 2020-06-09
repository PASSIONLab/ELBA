#!/bin/bash
git clone --branch seqan-v2.4.0 https://github.com/seqan/seqan.git
export SEQAN_HOME=$(pwd)/seqan/
cd $SEQAN_HOME
git checkout develop
