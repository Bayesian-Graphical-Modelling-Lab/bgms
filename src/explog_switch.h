#ifndef _EXPLOG_SWITCH_H_
#define	_EXPLOG_SWITCH_H_

/*
 *  To switch between these, try in R:
 *  Sys.setenv("PKG_CPPFLAGS" = "-DCUSTOM_EXP_LOG=0") # 0 for not using it, 1 for using it
 *  renv::install(".") # or whatever way you use to install bgms
 */

#if (defined(_WIN32) && (!defined(CUSTOM_EXP_LOG) || CUSTOM_EXP_LOG == 0)) \
  || (!defined(_WIN32) && defined(CUSTOM_EXP_LOG) && CUSTOM_EXP_LOG != 0)

#define USE_CUSTOM_LOG 1

#else

#define USE_CUSTOM_LOG 0

#endif

#if USE_CUSTOM_LOG

#include "e_exp.h"

#define MY_EXP __ieee754_exp
#define MY_LOG __ieee754_log

// TODO: add and use these
// #define MY_EXPM1 std::expm1
// #define MY_LOG1P std::log1p
// #define MY_LOG1PEXP


#else

#define MY_EXP std::exp
#define MY_LOG std::log

#endif

#endif