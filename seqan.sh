#!/bin/bash
wget https://github.com/seqan/seqan/releases/download/seqan-v2.4.0/seqan-library-2.4.0.zip
unzip seqan-library-2.4.0
rm seqan-library-2.4.0.zip
export SEQAN_HOME=seqan-library-2.4.0/

