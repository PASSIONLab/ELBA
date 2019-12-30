#!/bin/bash
cd build_release
cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_BUILD_TYPE=Release ../
make VERBOSE=1 -j4
