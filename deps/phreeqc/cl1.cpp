#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "Phreeqc.h"
#include "phqalloc.h"

/* debug
#define DEBUG_CL1
#define CHECK_ERRORS
 */

int Phreeqc::
cl1(int k, int l, int m, int n,
	int l_nklmd, int l_n2d,
	LDBLE * q,
	int *l_kode, LDBLE l_toler,
	int *l_iter, LDBLE * l_x, LDBLE * l_res, LDBLE * l_error,
	LDBLE * l_cu, int *l_iu, int *l_s, int check)
{
	/* System generated locals */
	union double_or_int
	{
		int ival;
		LDBLE dval;
	} *q2;

	/* Local variables */
	int nklm;
	LDBLE xmin, xmax;
	int iout = 0;
	// static i runs faster on windows
	register int i, j;
	register LDBLE l_z;
	int maxit, n1, n2;
	LDBLE pivot;
	int ia, ii, kk, nk, js;
	int in = 0;
	LDBLE sn;
	int iphase, kforce;
	LDBLE zu, zv;
	LDBLE tpivot;
	int klm, jmn, nkl, jpn;
	LDBLE cuv;
	long double sum;
	int klm1;
	int q_dim, cu_dim;
	int kode_arg;
	LDBLE check_toler;
#ifdef CHECK_ERRORS
	char **col_name, **row_name;
	int *row_back, *col_back;
#endif
/* THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX */
/* METHOD OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION */
/* TO A K BY N SYSTEM OF LINEAR EQUATIONS */
/*             AX=B */
/* SUBJECT TO L LINEAR EQUALITY CONSTRAINTS */
/*             CX=D */
/* AND M LINEAR INEQUALITY CONSTRAINTS */
/*             EX.LE.F. */
/* DESCRIPTION OF PARAMETERS */
/* K      NUMBER OF ROWS OF THE MATRIX A (K.GE.1). */
/* L      NUMBER OF ROWS OF THE MATRIX C (L.GE.0). */
/* M      NUMBER OF ROWS OF THE MATRIX E (M.GE.0). */
/* N      NUMBER OF COLUMNS OF THE MATRICES A,C,E (N.GE.1). */
/* KLMD   SET TO AT LEAST K+L+M FOR ADJUSTABLE DIMENSIONS. */
/* KLM2D  SET TO AT LEAST K+L+M+2 FOR ADJUSTABLE DIMENSIONS. */
/* NKLMD  SET TO AT LEAST N+K+L+M FOR ADJUSTABLE DIMENSIONS. */
/* N2D    SET TO AT LEAST N+2 FOR ADJUSTABLE DIMENSIONS */
/* Q      TWO DIMENSIONAL REAL ARRAY WITH KLM2D ROWS AND */
/*        AT LEAST N2D COLUMNS. */
/*        ON ENTRY THE MATRICES A,C AND E, AND THE VECTORS */
/*        B,D AND F MUST BE STORED IN THE FIRST K+L+M ROWS */
/*        AND N+1 COLUMNS OF Q AS FOLLOWS */
/*             A B */
/*         Q = C D */
/*             E F */
/*        THESE VALUES ARE DESTROYED BY THE SUBROUTINE. */
/* KODE   A CODE USED ON ENTRY TO, AND EXIT */
/*        FROM, THE SUBROUTINE. */
/*        ON ENTRY, THIS SHOULD NORMALLY BE SET TO 0. */
/*        HOWEVER, IF CERTAIN NONNEGATIVITY CONSTRAINTS */
/*        ARE TO BE INCLUDED IMPLICITLY, RATHER THAN */
/*        EXPLICITLY IN THE CONSTRAINTS EX.LE.F, THEN KODE */
/*        SHOULD BE SET TO 1, AND THE NONNEGATIVITY */
/*        CONSTRAINTS INCLUDED IN THE ARRAYS X AND */
/*        RES (SEE BELOW). */
/*        ON EXIT, KODE HAS ONE OF THE */
/*        FOLLOWING VALUES */
/*             0- OPTIMAL SOLUTION FOUND, */
/*             1- NO FEASIBLE SOLUTION TO THE */
/*                CONSTRAINTS, */
/*             2- CALCULATIONS TERMINATED */
/*                PREMATURELY DUE TO ROUNDING ERRORS, */
/*             3- MAXIMUM NUMBER OF ITERATIONS REACHED. */
/* TOLER  A SMALL POSITIVE TOLERANCE. EMPIRICAL */
/*        EVIDENCE SUGGESTS TOLER = 10**(-D*2/3), */
/*        WHERE D REPRESENTS THE NUMBER OF DECIMAL */
/*        DIGITS OF ACCURACY AVAILABLE. ESSENTIALLY, */
/*        THE SUBROUTINE CANNOT DISTINGUISH BETWEEN ZERO */
/*        AND ANY QUANTITY WHOSE MAGNITUDE DOES NOT EXCEED */
/*        TOLER. IN PARTICULAR, IT WILL NOT PIVOT ON ANY */
/*        NUMBER WHOSE MAGNITUDE DOES NOT EXCEED TOLER. */
/* ITER   ON ENTRY ITER MUST CONTAIN AN UPPER BOUND ON */
/*        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED. */
/*        A SUGGESTED VALUE IS 10*(K+L+M). ON EXIT ITER */
/*        GIVES THE NUMBER OF SIMPLEX ITERATIONS. */
/* X      ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST N2D. */
/*        ON EXIT THIS ARRAY CONTAINS A */
/*        SOLUTION TO THE L1 PROBLEM. IF KODE=1 */
/*        ON ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE */
/*        SIMPLE NONNEGATIVITY CONSTRAINTS ON THE */
/*        VARIABLES. THE VALUES -1, 0, OR 1 */
/*        FOR X(J) INDICATE THAT THE J-TH VARIABLE */
/*        IS RESTRICTED TO BE .LE.0, UNRESTRICTED, */
/*        OR .GE.0 RESPECTIVELY. */
/* RES    ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST KLMD. */
/*        ON EXIT THIS CONTAINS THE RESIDUALS B-AX */
/*        IN THE FIRST K COMPONENTS, D-CX IN THE */
/*        NEXT L COMPONENTS (THESE WILL BE =0),AND */
/*        F-EX IN THE NEXT M COMPONENTS. IF KODE=1 ON */
/*        ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE SIMPLE */
/*        NONNEGATIVITY CONSTRAINTS ON THE RESIDUALS */
/*        B-AX. THE VALUES -1, 0, OR 1 FOR RES(I) */
/*        INDICATE THAT THE I-TH RESIDUAL (1.LE.I.LE.K) IS */
/*        RESTRICTED TO BE .LE.0, UNRESTRICTED, OR .GE.0 */
/*        RESPECTIVELY. */
/* ERROR  ON EXIT, THIS GIVES THE MINIMUM SUM OF */
/*        ABSOLUTE VALUES OF THE RESIDUALS. */
/* CU     A TWO DIMENSIONAL REAL ARRAY WITH TWO ROWS AND */
/*        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE. */
/* IU     A TWO DIMENSIONAL INTEGER ARRAY WITH TWO ROWS AND */
/*        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE. */
/* S      INTEGER ARRAY OF SIZE AT LEAST KLMD, USED FOR */
/*        WORKSPACE. */
/*      DOUBLE PRECISION DBLE */
/*      REAL */

/* INITIALIZATION. */

	zv = 0;
	kode_arg = *l_kode;
	cl1_space(check, l_n2d, k + l + m, l_nklmd);

/* Parameter adjustments */
	q_dim = l_n2d;
	q2 = (union double_or_int *) q;
	cu_dim = l_nklmd;

/* Function Body */
	maxit = *l_iter;
	n1 = n + 1;
	n2 = n + 2;
	nk = n + k;
	nkl = nk + l;
	klm = k + l + m;
	klm1 = klm + 1;
	nklm = n + klm;
	kforce = 1;
	*l_iter = 0;
	js = 0;
	ia = -1;

/* SET UP LABELS IN Q. */
	for (j = 0; j < n; ++j)
	{
		q2[klm1 * q_dim + j].ival = j + 1;
	}
/* L10: */
	for (i = 0; i < klm; ++i)
	{
		q2[i * q_dim + n1].ival = n + i + 1;
		if (q2[i * q_dim + n].dval < 0.)
		{
			for (j = 0; j < n1; ++j)
			{
				q2[i * q_dim + j].dval = -q2[i * q_dim + j].dval;
			}
			q2[i * q_dim + n1].ival = -q2[i * q_dim + n1].ival;
/* L20: */
		}
	}
/* L30: */
/* SET UP PHASE 1 COSTS. */
	iphase = 2;
#ifdef DEBUG_CL1
	output_msg(sformatf( "Set up phase 1 costs\n"));
#endif
/* Zero first row of cu and iu */
	memcpy((void *) &(l_cu[0]), (void *) &(scratch[0]),
		   (size_t) nklm * sizeof(LDBLE));
	for (j = 0; j < nklm; ++j)
	{
		l_iu[j] = 0;
	}
/* L40: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L40\n"));
#endif
	if (l != 0)
	{
		for (j = nk; j < nkl; ++j)
		{
			l_cu[j] = 1.;
			l_iu[j] = 1;
		}
/* L50: */
		iphase = 1;
	}

/* Copy first row of cu and iu to second row */
	memcpy((void *) &(l_cu[cu_dim]), (void *) &(l_cu[0]),
		   (size_t) nklm * sizeof(LDBLE));
	memcpy((void *) &(l_iu[cu_dim]), (void *) &(l_iu[0]),
		   (size_t) nklm * sizeof(int));

/* L60: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L60\n"));
#endif
	if (m != 0)
	{
		for (j = nkl; j < nklm; ++j)
		{
			l_cu[cu_dim + j] = 1.;
			l_iu[cu_dim + j] = 1;
			jmn = j - n;
			if (q2[jmn * q_dim + n1].ival < 0)
			{
				iphase = 1;
			}
		}
/* L70: */
	}
/* L80: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L80\n"));
#endif
	if (*l_kode != 0)
	{
		for (j = 0; j < n; ++j)
		{
			if (l_x[j] < 0.)
			{
/* L90: */
				l_cu[j] = 1.;
				l_iu[j] = 1;
			}
			else if (l_x[j] > 0.)
			{
				l_cu[cu_dim + j] = 1.;
				l_iu[cu_dim + j] = 1;
			}
		}
/* L110: */
#ifdef DEBUG_CL1
		output_msg(sformatf( "L110\n"));
#endif
		for (j = 0; j < k; ++j)
		{
			jpn = j + n;
			if (l_res[j] < 0.)
			{
/* L120: */
				l_cu[jpn] = 1.;
				l_iu[jpn] = 1;
				if (q2[j * q_dim + n1].ival > 0)
				{
					iphase = 1;
				}
			}
			else if (l_res[j] > 0.)
			{
/* L130: */
				l_cu[cu_dim + jpn] = 1.;
				l_iu[cu_dim + jpn] = 1;
				if (q2[j * q_dim + n1].ival < 0)
				{
					iphase = 1;
				}
			}
		}
/* L140: */
	}
/* L150: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L150\n"));
#endif
	if (iphase == 2)
	{
		goto L500;
	}
/* COMPUTE THE MARGINAL COSTS. */
  L160:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L160\n"));
#endif
	for (j = js; j < n1; ++j)
	{
	sum = 0.;
		for (i = 0; i < klm; ++i)
		{
			ii = q2[i * q_dim + n1].ival;
			if (ii < 0)
			{
				l_z = l_cu[cu_dim - ii - 1];
			}
			else
			{
				l_z = l_cu[ii - 1];
			}
			sum += (long double) q2[i * q_dim + j].dval * (long double) l_z;
		}
		q2[klm * q_dim + j].dval = (double)sum;
	}
	for (j = js; j < n; ++j)
	{
		ii = q2[klm1 * q_dim + j].ival;
		if (ii < 0)
		{
			l_z = l_cu[cu_dim - ii - 1];
		}
		else
		{
			l_z = l_cu[ii - 1];
		}
		q2[klm * q_dim + j].dval -= l_z;
	}
/* DETERMINE THE VECTOR TO ENTER THE BASIS. */
  L240:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L240, xmax %e\n", xmax));
#endif
	xmax = 0.;
	if (js >= n)
	{
		goto L490;				/* test for optimality */
	}
	for (j = js; j < n; ++j)
	{
		zu = q2[klm * q_dim + j].dval;
		ii = q2[klm1 * q_dim + j].ival;
		if (ii > 0)
		{
			zv = -zu - l_cu[ii - 1] - l_cu[cu_dim + ii - 1];
		}
		else
		{
			ii = -ii;
			zv = zu;
			zu = -zu - l_cu[ii - 1] - l_cu[cu_dim + ii - 1];
		}
/* L260 */
		if (kforce == 1 && ii > n)
		{
			continue;
		}
		if (l_iu[ii - 1] != 1 && zu > xmax)
		{
			xmax = zu;
			in = j;
		}
/* L270 */
		if (l_iu[cu_dim + ii - 1] != 1 && zv > xmax)
		{
			xmax = zv;
			in = j;
		}
	}
/* L280 */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L280 xmax %e, toler %e\n", xmax, l_toler));
#endif
	if (xmax <= l_toler)
	{
#ifdef DEBUG_CL1
		output_msg(sformatf( "xmax before optimality test %e\n", xmax));
#endif
		goto L490;				/* test for optimality */
	}
	if (q2[klm * q_dim + in].dval != xmax)
	{
		for (i = 0; i < klm1; ++i)
		{
			q2[i * q_dim + in].dval = -q2[i * q_dim + in].dval;
		}
		q2[klm1 * q_dim + in].ival = -q2[klm1 * q_dim + in].ival;
/* L290: */
		q2[klm * q_dim + in].dval = xmax;
	}
/* DETERMINE THE VECTOR TO LEAVE THE BASIS. */
	if (iphase != 1 && ia != -1)
	{
		xmax = 0.;
/* find maximum absolute value in column "in" */
		for (i = 0; i <= ia; ++i)
		{
			l_z = fabs(q2[i * q_dim + in].dval);
			if (l_z > xmax)
			{
				xmax = l_z;
				iout = i;
			}
		}
/* L310: */
#ifdef DEBUG_CL1
		output_msg(sformatf( "L310, xmax %e\n", xmax));
#endif
/* switch row ia with row iout, use memcpy */
		if (xmax > l_toler)
		{
			if (ia != iout)
			{
				memcpy((void *) &(scratch[0]), (void *) &(q2[ia * q_dim]),
					   (size_t) n2 * sizeof(LDBLE));
				memcpy((void *) &(q2[ia * q_dim]), (void *) &(q2[iout * q_dim]),
					   (size_t) n2 * sizeof(LDBLE));
				memcpy((void *) &(q2[iout * q_dim]), (void *) &(scratch[0]),
					   (size_t) n2 * sizeof(LDBLE));
			}
/* L320: */
/* set pivot to row ia, column in */
			iout = ia;
			--ia;
			pivot = q2[iout * q_dim + in].dval;
			goto L420;			/* Gauss Jordan */
		}
	}
/* L330: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L330, xmax %e\n", xmax));
#endif
	kk = -1;
/* divide column n1 by positive value in column "in" greater than toler */
	for (i = 0; i < klm; ++i)
	{
		l_z = q2[i * q_dim + in].dval;
		if (l_z > l_toler)
		{
			++kk;
			l_res[kk] = q2[i * q_dim + n].dval / l_z;
			l_s[kk] = i;
		}
	}
/* L340: */
#ifdef DEBUG_CL1
	if (kk < 0)
	{
		output_msg(sformatf( "kode = 2 in loop 340.\n"));
	}
#endif
  L350:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L350, xmax %e\n", xmax));
#endif
	if (kk < 0)
	{
/* no positive value found in L340 or bypass intermediate verticies */
		*l_kode = 2;
		goto L590;
	}
/* L360: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L360, xmax %e\n", xmax));
#endif
/* find minimum residual */
	xmin = l_res[0];
	iout = l_s[0];
	j = 0;
	if (kk != 0)
	{
		for (i = 1; i <= kk; ++i)
		{
			if (l_res[i] < xmin)
			{
				j = i;
				xmin = l_res[i];
				iout = l_s[i];
			}
		}
/* L370: */
/* put kk in position j */
		l_res[j] = l_res[kk];
		l_s[j] = l_s[kk];
	}
/* L380: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L380 iout %d, xmin %e, xmax %e\n", iout,
			   xmin, xmax));
#endif
	--kk;
	pivot = q2[iout * q_dim + in].dval;
	ii = q2[iout * q_dim + n1].ival;
	if (iphase != 1)
	{
		if (ii < 0)
		{
/* L390: */
			if (l_iu[-ii - 1] == 1)
			{
				goto L420;
			}
		}
		else
		{
			if (l_iu[cu_dim + ii - 1] == 1)
			{
				goto L420;
			}
		}
	}
/* L400: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L400\n"));
#endif
	ii = abs(ii);
	cuv = l_cu[ii - 1] + l_cu[cu_dim + ii - 1];
	if (q2[klm * q_dim + in].dval - pivot * cuv > l_toler)
	{

/* BYPASS INTERMEDIATE VERTICES. */
		for (j = js; j < n1; ++j)
		{
			l_z = q2[iout * q_dim + j].dval;
			q2[klm * q_dim + j].dval -= l_z * cuv;
			q2[iout * q_dim + j].dval = -l_z;
		}
/* L410: */
		q2[iout * q_dim + n1].ival = -q2[iout * q_dim + n1].ival;
		goto L350;
	}
/* GAUSS-JORDAN ELIMINATION. */
  L420:
#ifdef DEBUG_CL1
	output_msg(sformatf( "Gauss Jordon %d\n", *l_iter));
#endif
	if (*l_iter >= maxit)
	{
		*l_kode = 3;
		goto L590;
	}
/* L430: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L430\n"));
#endif
	++(*l_iter);
	for (j = js; j < n1; ++j)
	{
		if (j != in)
		{
			q2[iout * q_dim + j].dval /= pivot;
		}
	}
/* L440: */
	for (j = js; j < n1; ++j)
	{
		if (j != in)
		{
			l_z = -q2[iout * q_dim + j].dval;
			for (i = 0; i < klm1; ++i)
			{
				if (i != iout)
				{
					q2[i * q_dim + j].dval += l_z * q2[i * q_dim + in].dval;
				}
			}
/* L450: */
		}
	}
/* L460: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L460\n"));
#endif
	tpivot = -pivot;
	for (i = 0; i < klm1; ++i)
	{
		if (i != iout)
		{
			q2[i * q_dim + in].dval /= tpivot;
		}
	}
/* L470: */
	q2[iout * q_dim + in].dval = 1. / pivot;
	ii = q2[iout * q_dim + n1].ival;
	q2[iout * q_dim + n1].ival = q2[klm1 * q_dim + in].ival;
	q2[klm1 * q_dim + in].ival = ii;
	ii = abs(ii);
	if (l_iu[ii - 1] == 0 || l_iu[cu_dim + ii - 1] == 0)
	{
		goto L240;
	}
/* switch column */
	for (i = 0; i < klm1; ++i)
	{
		l_z = q2[i * q_dim + in].dval;
		q2[i * q_dim + in].dval = q2[i * q_dim + js].dval;
		q2[i * q_dim + js].dval = l_z;
	}
	i = q2[klm1 * q_dim + in].ival;
	q2[klm1 * q_dim + in].ival = q2[klm1 * q_dim + js].ival;
	q2[klm1 * q_dim + js].ival = i;
/* L480: */
	++js;
	goto L240;
/* TEST FOR OPTIMALITY. */
  L490:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L490\n"));
#endif
	if (kforce == 0)
	{
		if (iphase == 1)
		{
			if (q2[klm * q_dim + n].dval <= l_toler)
			{
				goto L500;
			}
#ifdef DEBUG_CL1
			output_msg(sformatf( "q2[klm1-1, n1-1] > *toler. %e\n",
					   q2[(klm1 - 1) * q_dim + n1 - 1].dval));
#endif
			*l_kode = 1;
			goto L590;
		}
		*l_kode = 0;
		goto L590;
	}
	if (iphase != 1 || q2[klm * q_dim + n].dval > l_toler)
	{
		kforce = 0;
		goto L240;
	}
/* SET UP PHASE 2 COSTS. */
  L500:
#ifdef DEBUG_CL1
	output_msg(sformatf( "Set up phase 2 costs %d\n", *l_iter));
#endif
	iphase = 2;
	for (j = 0; j < nklm; ++j)
	{
		l_cu[j] = 0.;
	}
/* L510: */
	for (j = n; j < nk; ++j)
	{
		l_cu[j] = 1.;
	}
	memcpy((void *) &(l_cu[cu_dim]), (void *) &(l_cu[0]),
		   (size_t) nklm * sizeof(LDBLE));
/* L520: */
	for (i = 0; i < klm; ++i)
	{
		ii = q2[i * q_dim + n1].ival;
		if (ii <= 0)
		{
			if (l_iu[cu_dim - ii - 1] == 0)
			{
				continue;
			}
			l_cu[cu_dim - ii - 1] = 0.;
		}
		else
		{
/* L530: */
			if (l_iu[ii - 1] == 0)
			{
				continue;
			}
			l_cu[ii - 1] = 0.;
		}
/* L540: */
		++ia;
/* switch row */
		if (ia != i)
		{
			memcpy((void *) &(scratch[0]), (void *) &(q2[ia * q_dim]),
				   (size_t) n2 * sizeof(LDBLE));
			memcpy((void *) &(q2[ia * q_dim]), (void *) &(q2[i * q_dim]),
				   (size_t) n2 * sizeof(LDBLE));
			memcpy((void *) &(q2[i * q_dim]), (void *) &(scratch[0]),
				   (size_t) n2 * sizeof(LDBLE));
		}
/* L550: */
	}
/* L560: */
	goto L160;


/* PREPARE OUTPUT. */
  L590:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L590\n"));
#endif
	sum = 0.;
	for (j = 0; j < n; ++j)
	{
		l_x[j] = 0.;
	}
/* L600: */
	for (i = 0; i < klm; ++i)
	{
		l_res[i] = 0.;
	}
/* L610: */
	for (i = 0; i < klm; ++i)
	{
		ii = q2[i * q_dim + n1].ival;
		sn = 1.;
		if (ii < 0)
		{
			ii = -ii;
			sn = -1.;
		}
		if (ii <= n)
		{
/* L620: */
			l_x[ii - 1] = sn * q2[i * q_dim + n].dval;
		}
		else
		{
/* L630: */
			l_res[ii - n - 1] = sn * q2[i * q_dim + n].dval;
			if (ii >= n1 && ii <= nk)
			{
/*     *    DBLE(Q(I,N1)) */
			  sum += (long double) q2[i * q_dim + n].dval;
			}
		}
	}
/* L640: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L640\n"));
#endif
	*l_error = (double)sum;
	/*
	 *  Check calculation
	 */
	if ((check == 1) && (*l_kode == 0))
	{
		check_toler = 10. * l_toler;
		/*
		 *  Check optimization constraints
		 */
		if (kode_arg == 1)
		{
			for (i = 0; i < k; ++i)
			{
				if (res_arg[i] < 0.0)
				{
					if (l_res[i] > check_toler)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1: optimization constraint not satisfied row %d, res %s, constraint %f.\n",
								   row_name[row_back[i]], l_res[i], res_arg[i]));
#endif
						*l_kode = 1;
					}
				}
				else if (res_arg[i] > 0.0)
				{
					if (l_res[i] < -check_toler)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1: optimization constraint not satisfied row %s, res %e, constraint %f.\n",
								   row_name[row_back[i]], l_res[i], res_arg[i]));
#endif
						*l_kode = 1;
					}
				}
			}
		}
		/*
		 *  Check equalities
		 */
		for (i = k; i < k + l; ++i)
		{
			if (fabs(l_res[i]) > check_toler)
			{
#ifdef CHECK_ERRORS
				output_msg(sformatf(
						   "\tCL1: equality constraint not satisfied row %s, res %e, tolerance %e.\n",
						   row_name[row_back[i]], l_res[i], check_toler));
#endif
				*l_kode = 1;
			}
		}
		/*
		 *  Check inequalities
		 */
		for (i = k + l; i < k + l + m; ++i)
		{
			if (l_res[i] < -check_toler)
			{
#ifdef CHECK_ERRORS
				output_msg(sformatf(
						   "\tCL1: inequality constraint not satisfied row %s, res %e, tolerance %e.\n",
						   row_name[row_back[i]], l_res[i], check_toler));
#endif
				*l_kode = 1;
			}
		}
		/*
		 *   Check dissolution/precipitation constraints
		 */
		if (kode_arg == 1)
		{
			for (i = 0; i < n; ++i)
			{
				if (x_arg[i] < 0.0)
				{
					if (l_x[i] > check_toler)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1: dis/pre constraint not satisfied column %s, x %e, constraint %f.\n",
								   col_name[col_back[i]], l_x[i], x_arg[i]));
#endif
						*l_kode = 1;
					}
				}
				else if (x_arg[i] > 0.0)
				{
					if (l_x[i] < -check_toler)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1: dis/pre constraint not satisfied column %s, x %e, constraint %f.\n",
								   col_name[col_back[i]], l_x[i], x_arg[i]));
#endif
						*l_kode = 1;
					}
				}
			}
		}
		if (*l_kode == 1)
		{
			output_msg(sformatf(
					   "\n\tCL1: Roundoff errors in optimization.\n\t     Try using -multiple_precision in INVERSE_MODELING\n"));
		}
	}
	return 0;
}

void Phreeqc::
cl1_space(int check, int l_n2d, int klm, int l_nklmd)
{
	if (check == 1)
	{
		if (x_arg == NULL)
		{
			x_arg = (LDBLE *) PHRQ_malloc((size_t) (l_n2d * sizeof(LDBLE)));
		}
		else if (l_n2d > x_arg_max)
		{
			x_arg =
				(LDBLE *) PHRQ_realloc(x_arg, (size_t) (l_n2d * sizeof(LDBLE)));
			x_arg_max = l_n2d;
		}
		if (x_arg == NULL)
			malloc_error();
		zero_double(x_arg, l_n2d);

		if (res_arg == NULL)
		{
			res_arg = (LDBLE *) PHRQ_malloc((size_t) ((klm) * sizeof(LDBLE)));
		}
		else if (klm > res_arg_max)
		{
			res_arg =
				(LDBLE *) PHRQ_realloc(res_arg,
									   (size_t) ((klm) * sizeof(LDBLE)));
			res_arg_max = klm;
		}
		if (res_arg == NULL)
			malloc_error();
		zero_double(res_arg, klm);
	}

/* Make scratch space */
	if (scratch == NULL)
	{
		scratch = (LDBLE *) PHRQ_malloc((size_t) l_nklmd * sizeof(LDBLE));
	}
	else if (l_nklmd > scratch_max)
	{
		scratch =
			(LDBLE *) PHRQ_realloc(scratch, (size_t) l_nklmd * sizeof(LDBLE));
		scratch_max = l_nklmd;
	}
	if (scratch == NULL)
		malloc_error();
	zero_double(scratch, l_nklmd);
}
