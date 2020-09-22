/**************************************************************************
 *                                                                        *
 * File          : nvector.c                                              *
 * Programmers   : Radu Serban, LLNL                                      *
 * Version of    : 26 June 2002                                           *
 *------------------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California        *
 * Produced at the Lawrence Livermore National Laboratory                 *
 * All rights reserved                                                    *
 * For details, see LICENSE below                                         *
 *------------------------------------------------------------------------*
 * This is the implementation file for a generic NVECTOR                  *
 * package. It contains the implementation of the N_Vector                *
 * kernels listed in nvector.h.                                           *
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

#include "nvector.h"			/* generic M_Env and N_Vector */

N_Vector
N_VNew(integertype n, M_Env machEnv)
{
	N_Vector v_new;
	v_new = machEnv->ops->nvnew(n, machEnv);
	return (v_new);
}

N_Vector_S
N_VNew_S(integertype ns, integertype n, M_Env machEnv)
{
	N_Vector_S vs_new;
	vs_new = machEnv->ops->nvnewS(ns, n, machEnv);
	return (vs_new);
}

void
N_VFree(N_Vector v)
{
	v->menv->ops->nvfree(v);
}

void
N_VFree_S(integertype ns, N_Vector_S vs)
{
	(*vs)->menv->ops->nvfreeS(ns, vs);
}

N_Vector
N_VMake(integertype n, realtype * v_data, M_Env machEnv)
{
	N_Vector v_new;
	v_new = machEnv->ops->nvmake(n, v_data, machEnv);
	return (v_new);
}

void
N_VDispose(N_Vector v)
{
	v->menv->ops->nvdispose(v);
}

realtype *
N_VGetData(N_Vector v)
{
	realtype *data;
	data = v->menv->ops->nvgetdata(v);
	return (data);
}

void
N_VSetData(realtype * v_data, N_Vector v)
{
	v->menv->ops->nvsetdata(v_data, v);
}

void
N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
	z->menv->ops->nvlinearsum(a, x, b, y, z);
}

void
N_VConst(realtype c, N_Vector z)
{
	z->menv->ops->nvconst(c, z);
}

void
N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
	z->menv->ops->nvprod(x, y, z);
}

void
N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
	z->menv->ops->nvdiv(x, y, z);
}

void
N_VScale(realtype c, N_Vector x, N_Vector z)
{
	z->menv->ops->nvscale(c, x, z);
}

void
N_VAbs(N_Vector x, N_Vector z)
{
	z->menv->ops->nvabs(x, z);
}

void
N_VInv(N_Vector x, N_Vector z)
{
	z->menv->ops->nvinv(x, z);
}

void
N_VAddConst(N_Vector x, realtype b, N_Vector z)
{
	z->menv->ops->nvaddconst(x, b, z);
}

realtype
N_VDotProd(N_Vector x, N_Vector y)
{
	realtype prod;
	prod = y->menv->ops->nvdotprod(x, y);
	return (prod);
}

realtype
N_VMaxNorm(N_Vector x)
{
	realtype norm;
	norm = x->menv->ops->nvmaxnorm(x);
	return (norm);
}

realtype
N_VWrmsNorm(N_Vector x, N_Vector w)
{
	realtype norm;
	norm = x->menv->ops->nvwrmsnorm(x, w);
	return (norm);
}

realtype
N_VMin(N_Vector x)
{
	realtype minval;
	minval = x->menv->ops->nvmin(x);
	return (minval);
}

realtype
N_VWL2Norm(N_Vector x, N_Vector w)
{
	realtype norm;
	norm = x->menv->ops->nvwl2norm(x, w);
	return (norm);
}

realtype
N_VL1Norm(N_Vector x)
{
	realtype norm;
	norm = x->menv->ops->nvl1norm(x);
	return (norm);
}

void
N_VOneMask(N_Vector x)
{
	x->menv->ops->nvonemask(x);
}

void
N_VCompare(realtype c, N_Vector x, N_Vector z)
{
	z->menv->ops->nvcompare(c, x, z);
}

booleantype
N_VInvTest(N_Vector x, N_Vector z)
{
	booleantype flag;
	flag = z->menv->ops->nvinvtest(x, z);
	return (flag);
}

booleantype
N_VConstrProdPos(N_Vector c, N_Vector x)
{
	booleantype flag;
	flag = x->menv->ops->nvconstrprodpos(c, x);
	return (flag);
}

booleantype
N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
	booleantype flag;
	flag = x->menv->ops->nvconstrmask(c, x, m);
	return (flag);
}

realtype
N_VMinQuotient(N_Vector num, N_Vector denom)
{
	realtype quotient;
	quotient = num->menv->ops->nvminquotient(num, denom);
	return (quotient);
}

void
N_VPrint(N_Vector x)
{
	x->menv->ops->nvprint(x);
}
