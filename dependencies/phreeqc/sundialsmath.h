/**************************************************************************
 *                                                                        *
 * File          : sundialsmath.h                                         *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL            *
 * Version of    : 26 June 2002                                           *
 *------------------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California        *
 * Produced at the Lawrence Livermore National Laboratory                 *
 * All rights reserved                                                    *
 * For details, see LICENSE below                                         *
 *------------------------------------------------------------------------*
 * This is the header file for a C math library. The routines             *
 * listed here work with the type realtype as defined in                  *
 * sundialstypes.h.                                                       *
 * To do single precision floating point arithmetic, set the type         *
 * realtype to be float. To do double precision arithmetic, set           *
 * the type realtype to be double. The default implementations            *
 * for RPowerR and RSqrt call standard math library functions             *
 * which do double precision arithmetic. If this is unacceptable          *
 * when realtype is float, then the user should re-implement              *
 * these two routines by calling single precision routines                *
 * available on his/her machine.                                          *
 *                                                                        *
 *------------------------------------------------------------------------*
 * LICENSE                                                                *
 *------------------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California.       *
 * Produced at the Lawrence Livermore National Laboratory.                *
 * Written by S.D. Cohen, A.C. Hindmarsh, R. Serban,                      *
 *            D. Shumaker, and A.G. Taylor.                               *
 * UCRL-CODE-155951    (CVODE)                                            *
 * UCRL-CODE-155950    (CVODES)                                           *
 * UCRL-CODE-155952    (IDA)                                              *
 * UCRL-CODE-237203    (IDAS)                                             *
 * UCRL-CODE-155953    (KINSOL)                                           *
 * All rights reserved.                                                   *
 *                                                                        *
 * This file is part of SUNDIALS.                                         *
 *                                                                        *
 * Redistribution and use in source and binary forms, with or without     *
 * modification, are permitted provided that the following conditions     *
 * are met:                                                               *
 *                                                                        *
 * 1. Redistributions of source code must retain the above copyright      *
 * notice, this list of conditions and the disclaimer below.              *
 *                                                                        *
 * 2. Redistributions in binary form must reproduce the above copyright   *
 * notice, this list of conditions and the disclaimer (as noted below)    *
 * in the documentation and/or other materials provided with the          *
 * distribution.                                                          *
 *                                                                        *
 * 3. Neither the name of the UC/LLNL nor the names of its contributors   *
 * may be used to endorse or promote products derived from this software  *
 * without specific prior written permission.                             *
 *                                                                        *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      *
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS      *
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE         *
 * REGENTS OF THE UNIVERSITY OF CALIFORNIA, THE U.S. DEPARTMENT OF ENERGY *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,        *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT       *
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,  *
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY  *
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT    *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  *
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   *
 **************************************************************************/
#ifndef _sundialsmath_h
#define _sundialsmath_h

#include "sundialstypes.h"


/******************************************************************
 *                                                                *
 * Macros : MIN, MAX, ABS, SQR                                    *
 *----------------------------------------------------------------*
 * MIN(A, B) returns the minimum of A and B.                      *
 *                                                                *
 * MAX(A, B) returns the maximum of A and B.                      *
 *                                                                *
 * ABS(A) returns the absolute value of A.                        *
 *                                                                *
 * SQR(A) returns the square of A.                                *
 *                                                                *
 ******************************************************************/
#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef ABS
#define ABS(A)    ((A < 0) ? -(A) : (A))
#endif

#ifndef SQR
#define SQR(A)    ((A) * (A))
#endif

/******************************************************************
 *                                                                *
 * Function : UnitRoundoff                                        *
 * Usage    : realtype uround;                                    *
 *            uround = UnitRoundoff();                            *
 *----------------------------------------------------------------*
 * UnitRoundoff returns the unit roundoff u for real floating     *
 * point arithmetic, where u is defined to be the smallest        *
 * positive real such that 1.0 + u != 1.0.                        *
 *                                                                *
 ******************************************************************/

	realtype UnitRoundoff(void);


/******************************************************************
 *                                                                *
 * Function : RPowerI                                             *
 * Usage    : int exponent;                                       *
 *            realtype base, ans;                                 *
 *            ans = RPowerI(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerI returns the value base^exponent, where base is a       *
 * realtype and exponent is an int.                               *
 *                                                                *
 ******************************************************************/

	realtype RPowerI(realtype base, int exponent);


/******************************************************************
 *                                                                *
 * Function : RPowerR                                             *
 * Usage    : realtype base, exponent, ans;                       *
 *            ans = RPowerR(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerR returns the value base^exponent, where both base and   *
 * exponent are realtype. If base < 0.0, then RPowerR returns 0.0 *
 *                                                                *
 ******************************************************************/

	realtype RPowerR(realtype base, realtype exponent);


/******************************************************************
 *                                                                *
 * Function : RSqrt                                               *
 * Usage    : realtype sqrt_x;                                    *
 *            sqrt_x = RSqrt(x);                                  *
 *----------------------------------------------------------------*
 * RSqrt(x) returns the square root of x. If x < 0.0, then RSqrt  *
 * returns 0.0.                                                   *
 *                                                                *
 ******************************************************************/

	realtype RSqrt(realtype x);


#endif
