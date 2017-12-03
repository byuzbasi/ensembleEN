#include "autoconfig.h"

#ifndef DEBUG
#   define ARMA_NO_DEBUG
#endif

#ifndef HAVE_OPENMP_CXX
#   define ARMA_DONT_USE_OPENMP
#endif
