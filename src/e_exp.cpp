
/* @(#)e_exp.c 1.6 04/04/22 */
/*
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */


// this is a stripped version of
// https://github.com/JuliaMath/openlibm/blob/c9c6fd6eade60ba38427c675a7323b0d2bae6521/src/e_exp.c
// specifically, the header files have been reduced to include only the minimially necessary things
// throughout the file EDIT indicates a modification relative to the original source code

// EDIT: omitted header
// #include "cdefs-compat.h"

//__FBSDID("$FreeBSD: src/lib/msun/src/e_exp.c,v 1.14 2011/10/21 06:26:38 das Exp $");

/* __ieee754_exp(x)
 * Returns the exponential of x.
 *
 * Method
 *   1. Argument reduction:
 *      Reduce x to an r so that |r| <= 0.5*ln2 ~ 0.34658.
 *	Given x, find r and integer k such that
 *
 *               x = k*ln2 + r,  |r| <= 0.5*ln2.
 *
 *      Here r will be represented as r = hi-lo for better
 *	accuracy.
 *
 *   2. Approximation of exp(r) by a special rational function on
 *	the interval [0,0.34658]:
 *	Write
 *	    R(r**2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r**4/360 + ...
 *      We use a special Remes algorithm on [0,0.34658] to generate
 * 	a polynomial of degree 5 to approximate R. The maximum error
 *	of this polynomial approximation is bounded by 2**-59. In
 *	other words,
 *	    R(z) ~ 2.0 + P1*z + P2*z**2 + P3*z**3 + P4*z**4 + P5*z**5
 *  	(where z=r*r, and the values of P1 to P5 are listed below)
 *	and
 *	    |                  5          |     -59
 *	    | 2.0+P1*z+...+P5*z   -  R(z) | <= 2
 *	    |                             |
 *	The computation of exp(r) thus becomes
 *                             2*r
 *		exp(r) = 1 + -------
 *		              R - r
 *                                 r*R1(r)
 *		       = 1 + r + ----------- (for better accuracy)
 *		                  2 - R1(r)
 *	where
 *			         2       4             10
 *		R1(r) = r - (P1*r  + P2*r  + ... + P5*r   ).
 *
 *   3. Scale back to obtain exp(x):
 *	From step 1, we have
 *	   exp(x) = 2^k * exp(r)
 *
 * Special cases:
 *	exp(INF) is INF, exp(NaN) is NaN;
 *	exp(-INF) is 0, and
 *	for finite argument, only exp(0)=1 is exact.
 *
 * Accuracy:
 *	according to an error analysis, the error is always less than
 *	1 ulp (unit in the last place).
 *
 * Misc. info.
 *	For IEEE double
 *	    if x >  7.09782712893383973096e+02 then exp(x) overflow
 *	    if x < -7.45133219101941108420e+02 then exp(x) underflow
 *
 * Constants:
 * The hexadecimal values are the intended ones for the following
 * constants. The decimal values may be used, provided that the
 * compiler will convert from decimal to binary accurately enough
 * to produce the hexadecimal values shown.
 */

// EDIT: omitted headers, the necessary definitions are pasted below
// #include <openlibm_math.h>
// #include "math_private.h"

// EDIT: fpmath.h
/* Definitions provided directly by GCC and Clang. */
#if !(defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && defined(__ORDER_BIG_ENDIAN__))

#if defined(__GLIBC__)

#include <features.h>
#include <endian.h>
#define __ORDER_LITTLE_ENDIAN__  __LITTLE_ENDIAN
#define __ORDER_BIG_ENDIAN__     __BIG_ENDIAN
#define __BYTE_ORDER__           __BYTE_ORDER

#elif defined(__APPLE__)

#include <machine/endian.h>
#define __ORDER_LITTLE_ENDIAN__  LITTLE_ENDIAN
#define __ORDER_BIG_ENDIAN__     BIG_ENDIAN
#define __BYTE_ORDER__           BYTE_ORDER

#elif defined(_WIN32)

#define __ORDER_LITTLE_ENDIAN__  1234
#define __ORDER_BIG_ENDIAN__     4321
#define __BYTE_ORDER__           __ORDER_LITTLE_ENDIAN__

#endif

#endif /* __BYTE_ORDER__, __ORDER_LITTLE_ENDIAN__ and __ORDER_BIG_ENDIAN__ */

#ifndef __FLOAT_WORD_ORDER__
#define __FLOAT_WORD_ORDER__     __BYTE_ORDER__
#endif


#include <float.h>
// #include <stdint.h> // EDIT: for int32_t, u_int32_t
#include <cstdint> // EDIT: for u_int64_t

//EDIT: types-compat.h
typedef uint8_t               u_int8_t;
typedef uint16_t              u_int16_t;
typedef uint32_t              u_int32_t;
typedef uint64_t              u_int64_t;


/*
 * A union which permits us to convert between a double and two 32 bit
 * ints.
 */

// the defines below are from "math_private.h"
#if __FLOAT_WORD_ORDER__ == __ORDER_BIG_ENDIAN__

typedef union
{
  double value;
  struct
  {
    u_int32_t msw;
    u_int32_t lsw;
  } parts;
  struct
  {
    u_int64_t w;
  } xparts;
} ieee_double_shape_type;

#endif

#if __FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__

typedef union
{
  double value;
  struct
  {
    u_int32_t lsw;
    u_int32_t msw;
  } parts;
  struct
  {
    u_int64_t w;
  } xparts;
} ieee_double_shape_type;

#endif

/* Get the more significant 32 bit int from a double.  */

#define GET_HIGH_WORD(i,d)					    \
do {								                       \
  ieee_double_shape_type gh_u;					\
  gh_u.value = (d);						          \
  (i) = gh_u.parts.msw;						      \
} while (0)

/* Get the less significant 32 bit int from a double.  */

#define GET_LOW_WORD(i,d)					     \
do {								                       \
  ieee_double_shape_type gl_u;					\
  gl_u.value = (d);						          \
  (i) = gl_u.parts.lsw;						      \
} while (0)

/* Set a double from two 32 bit ints.  */

#define INSERT_WORDS(d,ix0,ix1)					\
do {								                        \
  ieee_double_shape_type iw_u;					 \
  iw_u.parts.msw = (ix0);					      \
  iw_u.parts.lsw = (ix1);					      \
  (d) = iw_u.value;						           \
} while (0)

/* Set the more significant 32 bits of a double from an int.  */

#define SET_HIGH_WORD(d,v)					    \
do {								                       \
  ieee_double_shape_type sh_u;					\
  sh_u.value = (d);						          \
  sh_u.parts.msw = (v);						      \
  (d) = sh_u.value;						          \
} while (0)                                                    \

/* Get two 32 bit ints from a double.  */

#define EXTRACT_WORDS(ix0,ix1,d)				\
do {								                        \
  ieee_double_shape_type ew_u;					 \
  ew_u.value = (d);						           \
  (ix0) = ew_u.parts.msw;					      \
  (ix1) = ew_u.parts.lsw;					      \
} while (0)


// EDIT: e_exp.c
static const double
one	= 1.0,
halF[2]	= {0.5,-0.5,},
huge	= 1.0e+300,
o_threshold=  7.09782712893383973096e+02,  /* 0x40862E42, 0xFEFA39EF */
u_threshold= -7.45133219101941108420e+02,  /* 0xc0874910, 0xD52D3051 */
ln2HI[2]   ={ 6.93147180369123816490e-01,  /* 0x3fe62e42, 0xfee00000 */
	     -6.93147180369123816490e-01,},/* 0xbfe62e42, 0xfee00000 */
ln2LO[2]   ={ 1.90821492927058770002e-10,  /* 0x3dea39ef, 0x35793c76 */
	     -1.90821492927058770002e-10,},/* 0xbdea39ef, 0x35793c76 */
invln2 =  1.44269504088896338700e+00, /* 0x3ff71547, 0x652b82fe */
P1   =  1.66666666666666019037e-01, /* 0x3FC55555, 0x5555553E */
P2   = -2.77777777770155933842e-03, /* 0xBF66C16C, 0x16BEBD93 */
P3   =  6.61375632143793436117e-05, /* 0x3F11566A, 0xAF25DE2C */
P4   = -1.65339022054652515390e-06, /* 0xBEBBBD41, 0xC5D26BF1 */
P5   =  4.13813679705723846039e-08; /* 0x3E663769, 0x72BEA4D0 */

static volatile double
twom1000= 9.33263618503218878990e-302;     /* 2**-1000=0x01700000,0*/

double __ieee754_exp(double x)	/* default IEEE double exp */
{
	double y,hi=0.0,lo=0.0,c,t,twopk;
	int32_t k=0,xsb;
	u_int32_t hx;

	GET_HIGH_WORD(hx,x);
	xsb = (hx>>31)&1;		/* sign bit of x */
	hx &= 0x7fffffff;		/* high word of |x| */

    /* filter out non-finite argument */
	if(hx >= 0x40862E42) {			/* if |x|>=709.78... */
            if(hx>=0x7ff00000) {
	        u_int32_t lx;
		GET_LOW_WORD(lx,x);
		if(((hx&0xfffff)|lx)!=0)
		     return x+x; 		/* NaN */
		else return (xsb==0)? x:0.0;	/* exp(+-inf)={inf,0} */
	    }
	    if(x > o_threshold) return huge*huge; /* overflow */
	    if(x < u_threshold) return twom1000*twom1000; /* underflow */
	}

        /* this implementation gives 2.7182818284590455 for exp(1.0),
           which is well within the allowable error. however,
           2.718281828459045 is closer to the true value so we prefer that
           answer, given that 1.0 is such an important argument value. */
        if (x == 1.0)
            return 2.718281828459045235360;

    /* argument reduction */
	if(hx > 0x3fd62e42) {		/* if  |x| > 0.5 ln2 */
	    if(hx < 0x3FF0A2B2) {	/* and |x| < 1.5 ln2 */
		hi = x-ln2HI[xsb]; lo=ln2LO[xsb]; k = 1-xsb-xsb;
	    } else {
		k  = (int)(invln2*x+halF[xsb]);
		t  = k;
		hi = x - t*ln2HI[0];	/* t*ln2HI is exact here */
		lo = t*ln2LO[0];
	    }
	    //STRICT_ASSIGN(double, x, hi - lo);
	    x = hi - lo; // EDIT: this is equivalent to the line above on platforms we care about
	}
	else if(hx < 0x3e300000)  {	/* when |x|<2**-28 */
	    if(huge+x>one) return one+x;/* trigger inexact */
	}
	else k = 0;

    /* x is now in primary range */
	t  = x*x;
	if(k >= -1021)
	    INSERT_WORDS(twopk,0x3ff00000+(k<<20), 0);
	else
	    INSERT_WORDS(twopk,0x3ff00000+((k+1000)<<20), 0);
	c  = x - t*(P1+t*(P2+t*(P3+t*(P4+t*P5))));
	if(k==0) 	return one-((x*c)/(c-2.0)-x);
	else 		y = one-((lo-(x*c)/(2.0-c))-hi);
	if(k >= -1021) {
	    if (k==1024) return y*2.0*0x1p1023;
	    return y*twopk;
	} else {
	    return y*twopk*twom1000;
	}
}

// #if (LDBL_MANT_DIG == 53)
// openlibm_weak_reference(exp, expl);
// #endif

// EDIT: e_log.c
static const double
  ln2_hi  =  6.93147180369123816490e-01,	/* 3fe62e42 fee00000 */
    ln2_lo  =  1.90821492927058770002e-10,	/* 3dea39ef 35793c76 */
    two54   =  1.80143985094819840000e+16,  /* 43500000 00000000 */
    Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
    Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
    Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
    Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
    Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
    Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
    Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

    static const double zero   =  0.0;

double __ieee754_log(double x)
    {
      double hfsq,f,s,z,R,w,t1,t2,dk;
      int32_t k,hx,i,j;
      u_int32_t lx;

      EXTRACT_WORDS(hx,lx,x);

      k=0;
      if (hx < 0x00100000) {			/* x < 2**-1022  */
  if (((hx&0x7fffffff)|lx)==0)
    return -two54/zero;		/* log(+-0)=-inf */
  if (hx<0) return (x-x)/zero;	/* log(-#) = NaN */
  k -= 54; x *= two54; /* subnormal number, scale up x */
  GET_HIGH_WORD(hx,x);
      }
      if (hx >= 0x7ff00000) return x+x;
      k += (hx>>20)-1023;
      hx &= 0x000fffff;
      i = (hx+0x95f64)&0x100000;
      SET_HIGH_WORD(x,hx|(i^0x3ff00000));	/* normalize x or x/2 */
  k += (i>>20);
  f = x-1.0;
  if((0x000fffff&(2+hx))<3) {	/* -2**-20 <= f < 2**-20 */
  if(f==zero) {
    if(k==0) {
      return zero;
    } else {
      dk=(double)k;
      return dk*ln2_hi+dk*ln2_lo;
    }
  }
  R = f*f*(0.5-0.33333333333333333*f);
  if(k==0) return f-R; else {dk=(double)k;
    return dk*ln2_hi-((R-dk*ln2_lo)-f);}
  }
  s = f/(2.0+f);
  dk = (double)k;
  z = s*s;
  i = hx-0x6147a;
  w = z*z;
  j = 0x6b851-hx;
  t1= w*(Lg2+w*(Lg4+w*Lg6));
  t2= z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7)));
  i |= j;
  R = t2+t1;
  if(i>0) {
    hfsq=0.5*f*f;
    if(k==0) return f-(hfsq-s*(hfsq+R)); else
      return dk*ln2_hi-((hfsq-(s*(hfsq+R)+dk*ln2_lo))-f);
  } else {
    if(k==0) return f-s*(f-R); else
      return dk*ln2_hi-((s*(f-R)-dk*ln2_lo)-f);
  }
}


#include "Rcpp.h"
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ieee754_exp(Rcpp::NumericVector x) {
  Rcpp::NumericVector y(x.size());
  for (int i = 0; i < x.size(); i++) {
    y[i] = __ieee754_exp(x[i]);
  }
  return y;
}

#include "Rcpp.h"
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ieee754_log(Rcpp::NumericVector x) {
  Rcpp::NumericVector y(x.size());
  for (int i = 0; i < x.size(); i++) {
    y[i] = __ieee754_log(x[i]);
  }
  return y;
}