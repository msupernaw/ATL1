#!/bin/bash
export TBBROOT="/home/matthew/NetBeansProjects/ATL/AutoDiff/third_party/tbb42_20140601oss" #
tbb_bin="/home/matthew/NetBeansProjects/ATL/AutoDiff/third_party/tbb42_20140601oss/build/linux_intel64_gcc_cc4.9.2_libc2.19_kernel3.16.0_release" #
if [ -z "$CPATH" ]; then #
    export CPATH="${TBBROOT}/include" #
else #
    export CPATH="${TBBROOT}/include:$CPATH" #
fi #
if [ -z "$LIBRARY_PATH" ]; then #
    export LIBRARY_PATH="${tbb_bin}" #
else #
    export LIBRARY_PATH="${tbb_bin}:$LIBRARY_PATH" #
fi #
if [ -z "$LD_LIBRARY_PATH" ]; then #
    export LD_LIBRARY_PATH="${tbb_bin}" #
else #
    export LD_LIBRARY_PATH="${tbb_bin}:$LD_LIBRARY_PATH" #
fi #
 #
