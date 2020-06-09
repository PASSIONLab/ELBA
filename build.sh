#!/bin/bash
mkdir build_release
cd build_release
cmake ..
make -j4
cd ..
