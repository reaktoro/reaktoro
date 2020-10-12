/**************************************************************************
 *                                                                        *
 * File          : dense.c                                                *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and                 *
 *                 Radu Serban @ LLNL                                     *
 * Version of    : 26 June 2002                                           *
 *------------------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California        *
 * Produced at the Lawrence Livermore National Laboratory                 *
 * All rights reserved                                                    *
 * For details, see LICENSE below                                         *
 *------------------------------------------------------------------------*
 * This is the implementation file for a generic DENSE linear             *
 * solver package.                                                        *
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
#include <stdlib.h>
#include "sundialstypes.h"
#include "sundialsmath.h"
#include "dense.h"
#include "smalldense.h"

/* WARNING don`t include any headers below here */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)


/* Implementation */


DenseMat
DenseAllocMat(integertype N)
{
	DenseMat A;

	if (N <= 0)
		return (NULL);

	A = (DenseMat) malloc(sizeof *A);
	if (A == NULL)
		return (NULL);

	A->data = denalloc(N);
	if (A->data == NULL)
	{
		free(A);
		return (NULL);
	}

	A->size = N;

	return (A);
}


integertype *
DenseAllocPiv(integertype N)
{
	if (N <= 0)
		return (NULL);

	return ((integertype *) malloc(N * sizeof(integertype)));
}


integertype
DenseFactor(DenseMat A, integertype * p)
{
	return (gefa(A->data, A->size, p));
}


void
DenseBacksolve(DenseMat A, integertype * p, realtype * b)
{
	gesl(A->data, A->size, p, b);
}


void
DenseZero(DenseMat A)
{
	denzero(A->data, A->size);
}

void
DenseCopy(DenseMat A, DenseMat B)
{
	dencopy(A->data, B->data, A->size);
}

void
DenseScale(realtype c, DenseMat A)
{
	denscale(c, A->data, A->size);
}

void
DenseAddI(DenseMat A)
{
	denaddI(A->data, A->size);
}

void
DenseFreeMat(DenseMat A)
{
	denfree(A->data);
	free(A);
}

void
DenseFreePiv(integertype * p)
{
	free(p);
}

void
DensePrint(DenseMat A)
{
	denprint(A->data, A->size);
}
