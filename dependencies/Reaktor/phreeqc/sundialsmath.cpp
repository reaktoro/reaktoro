/**************************************************************************
 *                                                                        *
 * File          : sundialsmath.c                                         *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL            *
 * Version of    : 26 June 2002                                           *
 *------------------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California        *
 * Produced at the Lawrence Livermore National Laboratory                 *
 * All rights reserved                                                    *
 * For details, see LICENSE below                                         *
 *------------------------------------------------------------------------*
 * This is the implementation file for a C math library.                  *
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


#include <stdio.h>
#include <math.h>
#include "sundialsmath.h"
#include "sundialstypes.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)


realtype
UnitRoundoff(void)
{
	realtype u;
	volatile realtype one_plus_u;

	u = ONE;
	one_plus_u = ONE + u;
	while (one_plus_u != ONE)
	{
		u /= TWO;
		one_plus_u = ONE + u;
	}
	u *= TWO;

	return (u);
}


realtype
RPowerI(realtype base, int exponent)
{
	int i, expt;
	realtype prod;

	prod = ONE;
	expt = ABS(exponent);
	for (i = 1; i <= expt; i++)
		prod *= base;
	if (exponent < 0)
		prod = ONE / prod;
	return (prod);
}


realtype
RPowerR(realtype base, realtype exponent)
{

	if (base <= ZERO)
		return (ZERO);

	/*  return((realtype)pow((double)base,(double)exponent)); */
	return ((realtype) pow(base, exponent));
}


realtype
RSqrt(realtype x)
{
	if (x <= ZERO)
		return (ZERO);

	/*  return((realtype) sqrt((double) x)); */
	return ((realtype) sqrt(x));
}
