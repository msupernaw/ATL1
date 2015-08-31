#!/bin/csh
setenv TBBROOT "/Users/matthewsupernaw/NetBeansProjects/SS4.0/support/tbb42_20140601oss" #
setenv tbb_bin "/Users/matthewsupernaw/NetBeansProjects/SS4.0/support/tbb42_20140601oss/build/macos_intel64_clang_cc4.2.1_os10.7.5_debug" #
if (! $?CPATH) then #
    setenv CPATH "${TBBROOT}/include" #
else #
    setenv CPATH "${TBBROOT}/include:$CPATH" #
endif #
if (! $?LIBRARY_PATH) then #
    setenv LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LIBRARY_PATH "${tbb_bin}:$LIBRARY_PATH" #
endif #
if (! $?DYLD_LIBRARY_PATH) then #
    setenv DYLD_LIBRARY_PATH "${tbb_bin}" #
else #
    setenv DYLD_LIBRARY_PATH "${tbb_bin}:$DYLD_LIBRARY_PATH" #
endif #
 #
