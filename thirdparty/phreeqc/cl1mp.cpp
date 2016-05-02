#ifdef INVERSE_CL1MP
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <gmp.h>

#include "Phreeqc.h"
#include "phqalloc.h"
#include "phrqtype.h"

/* debug
#define DEBUG_CL1
#define CHECK_ERRORS
 */


int Phreeqc::
cl1mp(int k, int l, int m, int n,
	  int nklmd, int n2d,
	  LDBLE * q_arg,
	  int *kode_arg, LDBLE toler_arg,
	  int *iter, LDBLE * x_arg, LDBLE * res_arg, LDBLE * error_arg,
	  LDBLE * cu_arg, int *iu, int *s, int check, LDBLE censor_arg)
{
	mpf_set_default_prec(256);
	/* System generated locals */
	union double_or_int
	{
		int ival;
		mpf_t dval;
	} *q2;

	/* Local variables */
	int nklm;
	int iout = 0;
	// static i runs faster
	int i, j;
	int maxit, n1; //, n2;
	int ia, ii, kk, nk, js;
	int in = 0;
	int iphase, kforce;
	int klm, jmn, nkl, jpn;
	int klm1;
	int *kode;
	int q_dim, cu_dim;
	int iswitch;
	mpf_t *q;
	mpf_t *x;
	mpf_t *res;
	mpf_t error;
	mpf_t *cu;
	mpf_t dummy, dummy1, sum, z, zu, zv, xmax, minus_one, toler, check_toler;
	/*mpf_t *scratch; */
	mpf_t pivot, xmin, cuv, tpivot, sn;
	mpf_t zero;
	int censor;
	mpf_t censor_tol;
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
	/*
	 *  mp variables
	 */
	censor = 1;
	if (censor_arg == 0.0)
		censor = 0;
	mpf_set_default_prec(96);
	mpf_init(zero);
	mpf_init(dummy);
	mpf_init(dummy1);
	mpf_init_set_d(censor_tol, censor_arg);
	q = (mpf_t *)
		PHRQ_malloc((size_t)
					(max_row_count * max_column_count * sizeof(mpf_t)));
	if (q == NULL)
		malloc_error();
	for (i = 0; i < max_row_count * max_column_count; ++i)
	{
		mpf_init_set_d(q[i], q_arg[i]);
		if (censor == 1)
		{
			if (mpf_cmp(q[i], zero) != 0)
			{
				mpf_abs(dummy1, q[i]);
				if (mpf_cmp(dummy1, censor_tol) <= 0)
				{
					mpf_set_si(q[i], 0);
				}
			}
		}
	}
	x = (mpf_t *) PHRQ_malloc((size_t) (n2d * sizeof(mpf_t)));
	if (x == NULL)
		malloc_error();
	for (i = 0; i < n2d; ++i)
	{
		mpf_init_set_d(x[i], x_arg[i]);
	}
	res = (mpf_t *) PHRQ_malloc((size_t) ((k + l + m) * sizeof(mpf_t)));
	if (res == NULL)
		malloc_error();
	for (i = 0; i < k + l + m; ++i)
	{
		mpf_init_set_d(res[i], res_arg[i]);
	}
	cu = (mpf_t *) PHRQ_malloc((size_t) (2 * nklmd * sizeof(mpf_t)));
	if (cu == NULL)
		malloc_error();
	for (i = 0; i < 2 * nklmd; ++i)
	{
		mpf_init_set_d(cu[i], cu_arg[i]);
	}
	kode = (int *) PHRQ_malloc(sizeof(int));
	if (kode == NULL)
		malloc_error();
	*kode = *kode_arg;
	mpf_init(sum);
	mpf_init(error);
	mpf_init(z);
	mpf_init(zu);
	mpf_init(zv);
	mpf_init(xmax);
	mpf_init_set_si(minus_one, -1);
	mpf_init_set_d(toler, toler_arg);
	mpf_init_set_d(check_toler, toler_arg);
	mpf_init(pivot);
	mpf_init(xmin);
	mpf_init(cuv);
	mpf_init(tpivot);
	mpf_init(sn);
/* Parameter adjustments */
	q_dim = n2d;
	q2 = (union double_or_int *) q;
	cu_dim = nklmd;

/* Function Body */
	maxit = *iter;
	n1 = n + 1;
	//	n2 = n + 2;
	nk = n + k;
	nkl = nk + l;
	klm = k + l + m;
	klm1 = klm + 1;
	nklm = n + klm;
	kforce = 1;
	*iter = 0;
	js = 0;
	ia = -1;
/* Make scratch space */
/*
	scratch = (LDBLE *) PHRQ_malloc( (size_t) nklmd * sizeof(LDBLE));
	if (scratch == NULL) malloc_error();
	for (i=0; i < nklmd; i++) {
		scratch[i] = 0.0;
	}
*/
/*
	scratch = (mpf_t *) PHRQ_malloc( (size_t) nklmd * sizeof(mpf_t));
	if (scratch == NULL) malloc_error();
	for (i=0; i < nklmd; i++) {
		mpf_init(scratch[i]);
	}
*/
/* SET UP LABELS IN Q. */
	for (j = 0; j < n; ++j)
	{
		q2[klm1 * q_dim + j].ival = j + 1;
	}
/* L10: */
	for (i = 0; i < klm; ++i)
	{
		q2[i * q_dim + n1].ival = n + i + 1;
		if (mpf_cmp_d(q2[i * q_dim + n].dval, 0.0) < 0)
		{
			for (j = 0; j < n1; ++j)
			{
				/* q2[ i * q_dim + j ].dval = -q2[ i * q_dim + j ].dval; */
				mpf_neg(q2[i * q_dim + j].dval, q2[i * q_dim + j].dval);
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
	/*memcpy( (void *) &(cu[0]), (void *) &(scratch[0]), (size_t) nklm * sizeof(mpf_t) ); */
	for (j = 0; j < nklm; ++j)
	{
		mpf_set_si(cu[j], 0);
		iu[j] = 0;
	}
/* L40: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L40\n"));
#endif
	if (l != 0)
	{
		for (j = nk; j < nkl; ++j)
		{
			mpf_set_si(cu[j], 1);
			/*cu[ j ] = 1.; */
			iu[j] = 1;
		}
/* L50: */
		iphase = 1;
	}

/* Copy first row of cu and iu to second row */
	/*memcpy( (void *) &(cu[cu_dim]), (void *) &(cu[0]), (size_t) nklm * sizeof(mpf_t) ); */
	for (i = 0; i < nklm; ++i)
	{
		mpf_set(cu[cu_dim + i], cu[i]);
	}
	memcpy((void *) &(iu[cu_dim]), (void *) &(iu[0]),
		   (size_t) nklm * sizeof(int));
/* L60: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L60\n"));
#endif
	if (m != 0)
	{
		for (j = nkl; j < nklm; ++j)
		{
			/* cu[ cu_dim + j ] = 1.; */
			mpf_set_si(cu[cu_dim + j], 1);
			iu[cu_dim + j] = 1;
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
	if (*kode != 0)
	{
		for (j = 0; j < n; ++j)
		{
			/* if ( x[j] < 0.) { */
			if (mpf_cmp_si(x[j], 0) < 0)
			{
/* L90: */
				/* cu[ j ] = 1.; */
				mpf_set_si(cu[j], 1);
				iu[j] = 1;
				/* } else if (x[j] > 0.) { */
			}
			else if (mpf_cmp_si(x[j], 0) > 0)
			{
				/* cu[ cu_dim + j ] = 1.; */
				mpf_set_si(cu[cu_dim + j], 1);
				iu[cu_dim + j] = 1;
			}
		}
/* L110: */
#ifdef DEBUG_CL1
		output_msg(sformatf( "L110\n"));
#endif
		for (j = 0; j < k; ++j)
		{
			jpn = j + n;
			/* if (res[j] < 0.) { */
			if (mpf_cmp_si(res[j], 0) < 0)
			{
/* L120: */
				/* cu[ jpn ] = 1.; */
				mpf_set_si(cu[jpn], 1);
				iu[jpn] = 1;
				if (q2[j * q_dim + n1].ival > 0)
				{
					iphase = 1;
				}
				/* } else if (res[j] > 0.) { */
			}
			else if (mpf_cmp_si(res[j], 0) > 0)
			{
/* L130: */
				/* cu[ cu_dim + jpn ] = 1.; */
				mpf_set_si(cu[cu_dim + jpn], 1);
				iu[cu_dim + jpn] = 1;
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
		mpf_set_si(sum, 0);
		for (i = 0; i < klm; ++i)
		{
			ii = q2[i * q_dim + n1].ival;
			if (ii < 0)
			{
				/* z = cu[ cu_dim - ii - 1 ]; */
				mpf_set(z, cu[cu_dim - ii - 1]);
			}
			else
			{
				/*z = cu[ ii - 1 ]; */
				mpf_set(z, cu[ii - 1]);
			}
			/*sum += q2[ i * q_dim + j ].dval * z; */
			mpf_mul(dummy, q2[i * q_dim + j].dval, z);
			mpf_add(sum, sum, dummy);
		}
		/*q2[ klm * q_dim + j ].dval = sum; */
		mpf_set(q2[klm * q_dim + j].dval, sum);
	}
	for (j = js; j < n; ++j)
	{
		ii = q2[klm1 * q_dim + j].ival;
		if (ii < 0)
		{
			/*z = cu[ cu_dim - ii - 1 ]; */
			mpf_set(z, cu[cu_dim - ii - 1]);
		}
		else
		{
			/*z = cu[ ii - 1 ]; */
			mpf_set(z, cu[ii - 1]);
		}
		/*q2[ klm * q_dim + j ].dval -= z; */
		mpf_sub(q2[klm * q_dim + j].dval, q2[klm * q_dim + j].dval, z);
	}
/* DETERMINE THE VECTOR TO ENTER THE BASIS. */
  L240:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L240, xmax %e\n", mpf_get_d(xmax)));
#endif
	/*xmax = 0.; */
	mpf_set_si(xmax, 0);
	if (js >= n)
	{
		goto L490;				/* test for optimality */
	}
	for (j = js; j < n; ++j)
	{
		/*zu = q2[ klm * q_dim + j ].dval; */
		mpf_set(zu, q2[klm * q_dim + j].dval);
		ii = q2[klm1 * q_dim + j].ival;
		if (ii > 0)
		{
			/*zv = -zu - cu[ ii - 1 ] - cu[ cu_dim + ii - 1 ]; */
			mpf_mul(dummy, cu[cu_dim + ii - 1], minus_one);
			mpf_sub(dummy, dummy, cu[ii - 1]);
			mpf_sub(zv, dummy, zu);
		}
		else
		{
			ii = -ii;
			/* zv = zu; */
			mpf_set(zv, zu);
			/* zu = -zu - cu[ ii - 1 ] - cu[ cu_dim + ii - 1 ]; */
			mpf_mul(dummy, cu[cu_dim + ii - 1], minus_one);
			mpf_sub(dummy, dummy, cu[ii - 1]);
			mpf_sub(zu, dummy, zu);
		}
/* L260 */
		if (kforce == 1 && ii > n)
		{
			continue;
		}
		/*if (iu[ ii - 1 ] != 1 && zu > xmax){ */
		if ((iu[ii - 1] != 1) && (mpf_cmp(zu, xmax) > 0))
		{
			/*xmax = zu; */
			mpf_set(xmax, zu);
			in = j;
		}
/* L270 */
		/*if (iu[ cu_dim + ii - 1 ] != 1 && zv > xmax ) { */
		if ((iu[cu_dim + ii - 1] != 1) && (mpf_cmp(zv, xmax) > 0))
		{
			/*xmax = zv; */
			mpf_set(xmax, zv);
			in = j;
		}
	}
/* L280 */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L280 xmax %e, toler %e\n", mpf_get_d(xmax),
			   mpf_get_d(toler)));
#endif
	/*if (xmax <= toler) { */
	if (mpf_cmp(xmax, toler) <= 0)
	{
#ifdef DEBUG_CL1
		output_msg(sformatf( "xmax before optimality test %e\n",
				   mpf_get_d(xmax)));
#endif
		goto L490;				/* test for optimality */
	}
	/*if (q2[ klm * q_dim + in ].dval != xmax) { */
	if (mpf_cmp(q2[klm * q_dim + in].dval, xmax) != 0)
	{
		for (i = 0; i < klm1; ++i)
		{
			/*q2[ i * q_dim + in ].dval = -q2[ i * q_dim + in ].dval; */
			mpf_neg(q2[i * q_dim + in].dval, q2[i * q_dim + in].dval);
		}
		q2[klm1 * q_dim + in].ival = -q2[klm1 * q_dim + in].ival;
/* L290: */
		/*q2[ klm * q_dim + in ].dval = xmax; */
		mpf_set(q2[klm * q_dim + in].dval, xmax);
	}
/* DETERMINE THE VECTOR TO LEAVE THE BASIS. */
	if (iphase != 1 && ia != -1)
	{
		/*xmax = 0.; */
		mpf_set_si(xmax, 0);
/* find maximum absolute value in column "in" */
		for (i = 0; i <= ia; ++i)
		{
			/*z = fabs(q2[ i * q_dim + in ].dval); */
			mpf_abs(z, q2[i * q_dim + in].dval);
			/*if (z > xmax) { */
			if (mpf_cmp(z, xmax) > 0)
			{
				/*xmax = z; */
				mpf_set(xmax, z);
				iout = i;
			}
		}
/* L310: */
#ifdef DEBUG_CL1
		output_msg(sformatf( "L310, xmax %e\n", mpf_get_d(xmax)));
#endif
/* switch row ia with row iout, use memcpy */
		/*if (xmax > toler) { */
		if (mpf_cmp(xmax, toler) > 0)
		{
			/*
			   memcpy( (void *) &(scratch[0]), (void *) &(q2[ ia * q_dim]),
			   (size_t) n2 * sizeof(mpf_t) );
			   memcpy( (void *) &(q2[ ia * q_dim ]), (void *) &(q2[ iout * q_dim]),
			   (size_t) n2 * sizeof(mpf_t) );
			   memcpy( (void *) &(q2[ iout * q_dim ]), (void *) &(scratch[ 0 ]),
			   (size_t) n2 * sizeof(mpf_t) );
			 */
			for (i = 0; i < n1; ++i)
			{
				mpf_set(dummy, q2[ia * q_dim + i].dval);
				mpf_set(q2[ia * q_dim + i].dval, q2[iout * q_dim + i].dval);
				mpf_set(q2[iout * q_dim + i].dval, dummy);
			}
			j = q2[ia * q_dim + n1].ival;
			q2[ia * q_dim + n1].ival = q2[iout * q_dim + n1].ival;
			q2[iout * q_dim + n1].ival = j;

/* L320: */
/* set pivot to row ia, column in */
			iout = ia;
			--ia;
			/*pivot = q2[ iout * q_dim + in ].dval; */
			mpf_set(pivot, q2[iout * q_dim + in].dval);
			goto L420;			/* Gauss Jordan */
		}
	}
/* L330: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L330, xmax %e\n", mpf_get_d(xmax)));
#endif
	kk = -1;
/* divide column n1 by positive value in column "in" greater than toler */
	for (i = 0; i < klm; ++i)
	{
		/*z = q2[ i * q_dim + in ].dval; */
		mpf_set(z, q2[i * q_dim + in].dval);
		/*if (z > toler) { */
		if (mpf_cmp(z, toler) > 0)
		{
			++kk;
			/*res[kk] = q2[ i * q_dim + n ].dval / z; */
			mpf_div(res[kk], q2[i * q_dim + n].dval, z);
			s[kk] = i;
		}
	}
/* L340: */
	if (kk < 0)
	{
		output_msg(sformatf( "kode = 2 in loop 340.\n"));
	}
  L350:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L350, xmax %e\n", mpf_get_d(xmax)));
#endif
	if (kk < 0)
	{
/* no positive value found in L340 or bypass intermediate verticies */
		*kode = 2;
		goto L590;
	}
/* L360: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L360, xmax %e\n", mpf_get_d(xmax)));
#endif
/* find minimum residual */
	/*xmin = res[ 0 ]; */
	mpf_set(xmin, res[0]);
	iout = s[0];
	j = 0;
	if (kk != 0)
	{
		for (i = 1; i <= kk; ++i)
		{
			/*if (res[i] < xmin) { */
			if (mpf_cmp(res[i], xmin) < 0)
			{
				j = i;
				/*xmin = res[i]; */
				mpf_set(xmin, res[i]);
				iout = s[i];
			}
		}
/* L370: */
/* put kk in position j */
		/*res[j] = res[kk]; */
		mpf_set(res[j], res[kk]);
		s[j] = s[kk];
	}
/* L380: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L380 iout %d, xmin %e, xmax %e\n", iout,
			   mpf_get_d(xmin), mpf_get_d(xmax)));
#endif
	--kk;
	/*pivot = q2[ iout * q_dim + in ].dval; */
	mpf_set(pivot, q2[iout * q_dim + in].dval);
	ii = q2[iout * q_dim + n1].ival;
	if (iphase != 1)
	{
		if (ii < 0)
		{
/* L390: */
			if (iu[-ii - 1] == 1)
			{
				goto L420;
			}
		}
		else
		{
			if (iu[cu_dim + ii - 1] == 1)
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
	/*cuv = cu[ ii - 1 ] + cu[ cu_dim + ii - 1]; */
	mpf_add(cuv, cu[ii - 1], cu[cu_dim + ii - 1]);
	/*if (q2[ klm * q_dim + in ].dval - pivot * cuv > toler) { */
	mpf_mul(dummy, pivot, cuv);
	mpf_sub(dummy, q2[klm * q_dim + in].dval, dummy);
	if (mpf_cmp(dummy, toler) > 0)
	{
/* BYPASS INTERMEDIATE VERTICES. */
		for (j = js; j < n1; ++j)
		{
			/*z = q2[ iout * q_dim + j ].dval; */
			mpf_set(z, q2[iout * q_dim + j].dval);
			/*q2[ klm * q_dim + j ].dval -= z * cuv; */
			mpf_mul(dummy1, z, cuv);
			mpf_sub(q2[klm * q_dim + j].dval, q2[klm * q_dim + j].dval,
					dummy1);

			if (censor == 1)
			{
				if (mpf_cmp(q2[klm * q_dim + j].dval, zero) != 0)
				{
					mpf_abs(dummy1, q2[klm * q_dim + j].dval);
					if (mpf_cmp(dummy1, censor_tol) <= 0)
					{
						mpf_set_si(q2[klm * q_dim + j].dval, 0);
					}
				}
			}

			/*q2[ iout * q_dim + j ].dval = -z; */
			mpf_neg(q2[iout * q_dim + j].dval, z);
		}
/* L410: */
		q2[iout * q_dim + n1].ival = -q2[iout * q_dim + n1].ival;
		goto L350;
	}
/* GAUSS-JORDAN ELIMINATION. */
  L420:
#ifdef DEBUG_CL1
	output_msg(sformatf( "Gauss Jordon %d\n", *iter));
#endif
	if (*iter >= maxit)
	{
		*kode = 3;
		goto L590;
	}
/* L430: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L430\n"));
#endif
	++(*iter);
	for (j = js; j < n1; ++j)
	{
		if (j != in)
		{
			/*q2[ iout * q_dim + j ].dval /= pivot; */
			mpf_div(q2[iout * q_dim + j].dval, q2[iout * q_dim + j].dval,
					pivot);
		}
	}
/* L440: */
	for (j = js; j < n1; ++j)
	{
		if (j != in)
		{
			/*z = -q2[ iout * q_dim + j ].dval; */
			mpf_neg(z, q2[iout * q_dim + j].dval);
			for (i = 0; i < klm1; ++i)
			{
				if (i != iout)
				{
					/*q2[ i * q_dim + j ].dval += z * q2[ i * q_dim + in ].dval; */
					mpf_mul(dummy, z, q2[i * q_dim + in].dval);
					mpf_add(q2[i * q_dim + j].dval, q2[i * q_dim + j].dval,
							dummy);

					if (censor == 1)
					{
						if (mpf_cmp(q2[i * q_dim + j].dval, zero) != 0)
						{
							mpf_abs(dummy1, q2[i * q_dim + j].dval);
							if (mpf_cmp(dummy1, censor_tol) <= 0)
							{
								mpf_set_si(q2[i * q_dim + j].dval, 0);
							}
						}
					}
				}
			}
/* L450: */
		}
	}
/* L460: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L460\n"));
#endif
	/*tpivot = -pivot; */
	mpf_neg(tpivot, pivot);
	for (i = 0; i < klm1; ++i)
	{
		if (i != iout)
		{
			/*q2[ i * q_dim + in ].dval /= tpivot; */
			mpf_div(q2[i * q_dim + in].dval, q2[i * q_dim + in].dval, tpivot);
		}
	}
/* L470: */
	/*q2[ iout * q_dim + in ].dval = 1. / pivot; */
	mpf_set_si(dummy, 1);
	mpf_div(q2[iout * q_dim + in].dval, dummy, pivot);
	ii = q2[iout * q_dim + n1].ival;
	q2[iout * q_dim + n1].ival = q2[klm1 * q_dim + in].ival;
	q2[klm1 * q_dim + in].ival = ii;
	ii = abs(ii);
	if (iu[ii - 1] == 0 || iu[cu_dim + ii - 1] == 0)
	{
		goto L240;
	}
/* switch column */
	for (i = 0; i < klm1; ++i)
	{
		/*z = q2[ i * q_dim + in ].dval; */
		mpf_set(z, q2[i * q_dim + in].dval);
		/*q2[ i * q_dim + in ].dval = q2[ i * q_dim + js ].dval; */
		mpf_set(q2[i * q_dim + in].dval, q2[i * q_dim + js].dval);
		/*q2[ i * q_dim + js ].dval = z; */
		mpf_set(q2[i * q_dim + js].dval, z);
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
			/*if (q2[ klm * q_dim + n ].dval <= toler) { */
			if (mpf_cmp(q2[klm * q_dim + n].dval, toler) <= 0)
			{
				goto L500;
			}
#ifdef DEBUG_CL1
			output_msg(sformatf( "q2[klm1-1, n1-1] > *toler. %e\n",
					   mpf_get_d(q2[(klm1 - 1) * q_dim + n1 - 1].dval)));
#endif
			*kode = 1;
			goto L590;
		}
		*kode = 0;
		goto L590;
	}
	/*if (iphase != 1 || q2[ klm * q_dim + n ].dval > toler) { */
	if ((iphase != 1) || (mpf_cmp(q2[klm * q_dim + n].dval, toler) > 0))
	{
		kforce = 0;
		goto L240;
	}
/* SET UP PHASE 2 COSTS. */
  L500:
#ifdef DEBUG_CL1
	output_msg(sformatf( "Set up phase 2 costs %d\n", *iter));
#endif
	iphase = 2;
	for (j = 0; j < nklm; ++j)
	{
		/*cu[ j ] = 0.; */
		mpf_set_si(cu[j], 0);
	}
/* L510: */
	for (j = n; j < nk; ++j)
	{
		/*cu[ j ] = 1.; */
		mpf_set_si(cu[j], 1);
	}
	/*
	   memcpy( (void *) &(cu[cu_dim]), (void *) &(cu[0]), (size_t) nklm * sizeof(LDBLE) );
	 */
	for (i = 0; i < nklm; ++i)
	{
		mpf_set(cu[cu_dim + i], cu[i]);
	}

/* L520: */
	for (i = 0; i < klm; ++i)
	{
		ii = q2[i * q_dim + n1].ival;
		if (ii <= 0)
		{
			if (iu[cu_dim - ii - 1] == 0)
			{
				continue;
			}
			/*cu[ cu_dim - ii - 1 ] = 0.; */
			mpf_set_si(cu[cu_dim - ii - 1], 0);
		}
		else
		{
/* L530: */
			if (iu[ii - 1] == 0)
			{
				continue;
			}
			/*cu[ ii - 1 ] = 0.; */
			mpf_set_si(cu[ii - 1], 0);
		}
/* L540: */
		++ia;
/* switch row */
		/*
		   memcpy( (void *) &(scratch[0]), (void *) &(q2[ ia * q_dim]),
		   (size_t) n2 * sizeof(LDBLE) );
		   memcpy( (void *) &(q2[ ia * q_dim ]), (void *) &(q2[ i * q_dim]),
		   (size_t) n2 * sizeof(LDBLE) );
		   memcpy( (void *) &(q2[ i * q_dim ]), (void *) &(scratch[ 0 ]),
		   (size_t) n2 * sizeof(LDBLE) );
		 */
		for (iswitch = 0; iswitch < n1; ++iswitch)
		{
			mpf_set(dummy, q2[ia * q_dim + iswitch].dval);
			mpf_set(q2[ia * q_dim + iswitch].dval,
					q2[i * q_dim + iswitch].dval);
			mpf_set(q2[i * q_dim + iswitch].dval, dummy);
		}
		iswitch = q2[ia * q_dim + n1].ival;
		q2[ia * q_dim + n1].ival = q2[i * q_dim + n1].ival;
		q2[i * q_dim + n1].ival = iswitch;
/* L550: */
	}
/* L560: */
	goto L160;


/* PREPARE OUTPUT. */
  L590:
#ifdef DEBUG_CL1
	output_msg(sformatf( "L590\n"));
#endif
	/*sum = 0.; */
	mpf_set_si(sum, 0);
	for (j = 0; j < n; ++j)
	{
		/*x[j] = 0.; */
		mpf_set_si(x[j], 0);
	}
/* L600: */
	for (i = 0; i < klm; ++i)
	{
		/*res[i] = 0.; */
		mpf_set_si(res[i], 0);
	}
/* L610: */
	for (i = 0; i < klm; ++i)
	{
		ii = q2[i * q_dim + n1].ival;
		/*sn = 1.; */
		mpf_set_si(sn, 1);
		if (ii < 0)
		{
			ii = -ii;
			/*sn = -1.; */
			mpf_set_si(sn, -1);
		}
		if (ii <= n)
		{
/* L620: */
			/*x[ii - 1] = sn * q2[ i * q_dim + n ].dval; */
			mpf_mul(x[ii - 1], sn, q2[i * q_dim + n].dval);
		}
		else
		{
/* L630: */
			/*res[ii - n - 1] = sn * q2[ i * q_dim + n ].dval; */
			mpf_mul(res[ii - n - 1], sn, q2[i * q_dim + n].dval);
			if (ii >= n1 && ii <= nk)
			{
/*     *    DBLE(Q(I,N1)) */
				/*sum += q2[ i * q_dim + n ].dval; */
				mpf_add(sum, sum, q2[i * q_dim + n].dval);
			}
		}
	}
/* L640: */
#ifdef DEBUG_CL1
	output_msg(sformatf( "L640\n"));
#endif
	/*
	 *  Check calculation
	 */
	mpf_set_si(dummy, 100);
	mpf_mul(check_toler, toler, dummy);
	if (check && *kode == 0)
	{
		/*
		 *  Check optimization constraints
		 */
		if (*kode_arg == 1)
		{
			for (i = 0; i < k; ++i)
			{
				if (res_arg[i] < 0.0)
				{
					mpf_sub(dummy, res[i], check_toler);
					mpf_set_si(dummy1, 0);
					if (mpf_cmp(dummy, dummy1) > 0)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1MP: optimization constraint not satisfied row %d, res %e, constraint %f.\n",
								   i, mpf_get_d(res[i]), res_arg[i]));
#endif
						*kode = 1;
					}
				}
				else if (res_arg[i] > 0.0)
				{
					mpf_add(dummy, res[i], check_toler);
					mpf_set_si(dummy1, 0);
					if (mpf_cmp(dummy, dummy1) < 0)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1MP: optimization constraint not satisfied row %d, res %e, constraint %f.\n",
								   i, mpf_get_d(res[i]), res_arg[i]));
#endif
						*kode = 1;
					}
				}
			}
		}
		/*
		 *  Check equalities
		 */
		for (i = k; i < k + l; ++i)
		{
			mpf_abs(dummy, res[i]);
			if (mpf_cmp(dummy, check_toler) > 0)
			{
#ifdef CHECK_ERRORS
				output_msg(sformatf(
						   "\tCL1MP: equality constraint not satisfied row %d, res %e, tolerance %e.\n",
						   i, mpf_get_d(res[i]), mpf_get_d(check_toler)));
#endif

				*kode = 1;
			}
		}
		/*
		 *  Check inequalities
		 */
		for (i = k + l; i < k + l + m; ++i)
		{
			mpf_neg(dummy, check_toler);
			if (mpf_cmp(res[i], dummy) < 0)
			{
#ifdef CHECK_ERRORS
				output_msg(sformatf(
						   "\tCL1MP: inequality constraint not satisfied row %d, res %e, tolerance %e.\n",
						   i, mpf_get_d(res[i]), mpf_get_d(check_toler)));
#endif
				*kode = 1;
			}
		}
		/*
		 *   Check dissolution/precipitation constraints
		 */
		if (*kode_arg == 1)
		{
			for (i = 0; i < n; ++i)
			{
				if (x_arg[i] < 0.0)
				{
					mpf_sub(dummy, x[i], check_toler);
					mpf_set_si(dummy1, 0);
					if (mpf_cmp(dummy, dummy1) > 0)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1MP: dis/pre constraint not satisfied column %d, x %e, constraint %f.\n",
								   i, mpf_get_d(x[i]), x_arg[i]));
#endif
						*kode = 1;
					}
				}
				else if (x_arg[i] > 0.0)
				{
					mpf_add(dummy, x[i], check_toler);
					mpf_set_si(dummy1, 0);
					if (mpf_cmp(dummy, dummy1) < 0)
					{
#ifdef CHECK_ERRORS
						output_msg(sformatf(
								   "\tCL1MP: dis/pre constraint not satisfied column %d, x %e, constraint %f.\n",
								   i, mpf_get_d(x[i]), x_arg[i]));
#endif
						*kode = 1;
					}
				}
			}
		}
		if (*kode == 1)
		{
			output_msg(sformatf(
					   "\n\tCL1MP: Roundoff errors in optimization.\n\t       Deleting model.\n"));
		}
	}
	/*
	 * set return variables
	 */
	/**error = sum;*/
	mpf_set(error, sum);
	*error_arg = mpf_get_d(error);
	*kode_arg = *kode;
	for (i = 0; i < n2d; ++i)
	{
		x_arg[i] = mpf_get_d(x[i]);
	}
	for (i = 0; i < k + l + m; ++i)
	{
		res_arg[i] = mpf_get_d(res[i]);
	}

	/*scratch = free_check_null (scratch); */

	for (i = 0; i < max_row_count * max_column_count; ++i)
	{
		mpf_clear(q[i]);
	}
	q = (mpf_t *) free_check_null(q);
	for (i = 0; i < n2d; ++i)
	{
		mpf_clear(x[i]);
	}
	x = (mpf_t *) free_check_null(x);
	for (i = 0; i < k + l + m; ++i)
	{
		mpf_clear(res[i]);
	}
	res = (mpf_t *) free_check_null(res);
	for (i = 0; i < 2 * nklmd; ++i)
	{
		mpf_clear(cu[i]);
	}
	cu = (mpf_t *) free_check_null(cu);
	mpf_clear(zero);
	mpf_clear(dummy);
	mpf_clear(dummy1);
	mpf_clear(sum);
	mpf_clear(error);
	mpf_clear(z);
	mpf_clear(zu);
	mpf_clear(zv);
	mpf_clear(xmax);
	mpf_clear(minus_one);
	mpf_clear(toler);
	mpf_clear(check_toler);
	mpf_clear(pivot);
	mpf_clear(xmin);
	mpf_clear(cuv);
	mpf_clear(tpivot);
	mpf_clear(sn);
	mpf_clear(censor_tol);
	kode = (int *) free_check_null(kode);
	return 0;
}
#endif // INVERSE_CL1MP
