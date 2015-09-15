#ifndef ATL_PLATFORM
#define ATL_PLATFORM


#if defined(linux) || defined(__linux) || defined(__linux__)
#  define ATL_LINUX
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__DragonFly__)
#  define ATL_BSD
#elif defined(sun) || defined(__sun)
#  define ATL_SOLARIS
#elif defined(__sgi)
#  define ATL_IRIX
#elif defined(__hpux)
#  define ATL_HPUX
#elif defined(__CYGWIN__)
#  define ATL_CYGWIN
#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#  define ATL_WIN32
#elif defined(_WIN64) || defined(__WIN64__) || defined(WIN64)
#  define ATL_WIN64
#elif defined(__BEOS__)
#  define ATL_BEOS
#elif defined(macintosh) || defined(__APPLE__) || defined(__APPLE_CC__)
#  define ATL_MACOS
#elif defined(__IBMCPP__) || defined(_AIX)
#  define ATL_AIX
#elif defined(__amigaos__)
#  define ATL_AMIGAOS
#elif defined(__QNXNTO__)
#  define ATL_QNXNTO
#endif

#endif