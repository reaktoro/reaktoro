#!/bin/bash
set -ev

conda env create
set +v
source activate reaktoro

# Number of compilation jobs is limited by machine RAM. Too many and GCC will die with
# "internal compiler error"
export "NUMBER_OF_COMPILATION_JOBS=3"

export "CFLAGS=-I$CONDA_PREFIX/include $CFLAGS"
export "CXXFLAGS=-I$CONDA_PREFIX/include $CXXFLAGS"

echo -e "\n\nCurrent compiler configuration:\nCXX=$CXX\nCC=$CC\nCFLAGS=$CFLAGS\nCXXFLAGS=$CXXFLAGS"

echo -e "\nCurrent packages:"
set -v
conda list

mkdir -p build
mkdir -p artifacts
cd build

cmake -G Ninja -DBUILD_ALL=ON -DCMAKE_BUILD_TYPE=Release -DPYTHON_DIR=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=../artifacts -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON "-DTHIRDPARTY_COMMON_ARGS=-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON" "-DTHIRDPARTY_EXTRA_BUILD_ARGS=-- -j $NUMBER_OF_COMPILATION_JOBS" ..

cmake --build . --target install -- -j $NUMBER_OF_COMPILATION_JOBS

ls -la
ls -la ../artifacts
