#!/bin/csh
setenv TBBROOT "/home/matthew/NetBeansProjects/ATL/AutoDiff/third_party/tbb42_20140601oss" #
setenv tbb_bin "/home/matthew/NetBeansProjects/ATL/AutoDiff/third_party/tbb42_20140601oss/build/linux_intel64_gcc_cc4.9.2_libc2.19_kernel3.16.0_release" #
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
if (! $?LD_LIBRARY_PATH) then #
    setenv LD_LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LD_LIBRARY_PATH "${tbb_bin}:$LD_LIBRARY_PATH" #
endif #
 #
