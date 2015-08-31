#!/bin/bash
export TBBROOT="/Users/matthewsupernaw/NetBeansProjects/SS4.0/support/tbb42_20140601oss" #
tbb_bin="/Users/matthewsupernaw/NetBeansProjects/SS4.0/support/tbb42_20140601oss/build/macos_intel64_clang_cc4.2.1_os10.7.5_debug" #
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
if [ -z "$DYLD_LIBRARY_PATH" ]; then #
    export DYLD_LIBRARY_PATH="${tbb_bin}" #
else #
    export DYLD_LIBRARY_PATH="${tbb_bin}:$DYLD_LIBRARY_PATH" #
fi #
 #
