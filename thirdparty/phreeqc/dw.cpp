#include "Phreeqc.h"


/* ---------------------------------------------------------------------- */
int Phreeqc::
DW(LDBLE T)
/* ---------------------------------------------------------------------- */
/*
C
C      SUBROUTINE TO CALCULATE THE DENSITY OF WATER AS A FUNCTION OF
C      TEMPERATURE.  T IS IN KELVIN, P IS IN PASCALS, DW0 IS IN G/CM^3
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	LDBLE FP = 9.869232667e0, P, DGSS, D;

	BB(T);
	P = 1.0e0 / FP;
	if (T > 373.149e0)
		P = PS(T);
	DGSS = P / T / .4e0;
	if (T < TZ)
	{
		DGSS = 1.0e0 / (VLEST(T));
	}
	DFIND(&D, P, DGSS, T);
	DW0 = D;
	VP = P * FP;
	return OK;
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
BB(LDBLE T)
/* ---------------------------------------------------------------------- */
/*

C
C     THIS SUBROUTINE CALCULATES THE B'S NEEDED FOR FUNCTION DW.
C     THE B'S CALCULATED HERE ARE IN CM3/G.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	LDBLE V[11];
	int I;
	/* COMMON /BCONST/ */
	LDBLE P[11] =
		{ 0, 0.7478629e0, -.3540782e0, 0.e0, 0e0, .007159876e0, 0.e0,
		-.003528426e0, 0., 0., 0.
	};
	LDBLE Q[11] =
		{ 0, 1.1278334e0, 0.e0, -.5944001e0, -5.010996e0, 0.e0, .63684256e0,
		0., 0., 0., 0.
	};

	V[1] = 1.0;
	for (I = 2; I <= 10; I++)
	{
		V[I] = V[I - 1] * TZ / T;
	}
	B1 = P[1] + P[2] * log(1.e0 / V[2]);
	B2 = Q[1];
	B1T = P[2] * V[2] / TZ;
	B2T = 0.e0;
	B1TT = 0.e0;
	B2TT = 0.e0;
	for (I = 3; I <= 10; I++)
	{
		B1 = B1 + P[I] * V[I - 1];
		B2 = B2 + Q[I] * V[I - 1];
		B1T = B1T - (I - 2) * P[I] * V[I - 1] / T;
		B2T = B2T - (I - 2) * Q[I] * V[I - 1] / T;
		B1TT = B1TT + P[I] * (I - 2) * (I - 2) * V[I - 1] / T / T;
		B2TT = B2TT + Q[I] * (I - 2) * (I - 2) * V[I - 1] / T / T;
	}
	B1TT = B1TT - B1T / T;
	B2TT = B2TT - B2T / T;
	return OK;
}

/* ---------------------------------------------------------------------- */
 LDBLE Phreeqc::
PS(LDBLE T)
/* ---------------------------------------------------------------------- */
/*
C
C     THIS FUNCTION CALCULATES AN APPROXIMATION TO THE VAPOR PRESSURE, P
C     AS A FUNCTION OF THE INPUT TEMPERATURE. THE VAPOR PRESSURE
C     CALCULATED AGREES WITH THE VAPOR PRESSURE PREDICTED BY THE SURFACE
C     TO WITHIN .02% TO WITHIN A DEGREE OR SO OF THE CRITICAL TEMPERATUR
C     AND CAN SERVE AS AN INITIAL GUESS FOR FURTHER REFINEMENT BY
C     IMPOSING THE CONDITION THAT GL=GV.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	LDBLE A[9] = { 0, -7.8889166e0, 2.5514255e0, -6.716169e0,
		33.239495e0, -105.38479e0, 174.35319e0, -148.39348e0,
		48.631602e0
	};
	LDBLE PL, V, W, B, L_Z, Q;
	int I;
	if (T <= 314.e0)
	{
		PL = 6.3573118e0 - 8858.843e0 / T + 607.56335e0 * pow(T, (LDBLE) -.6e0);
		return (.1e0 * exp(PL));
	}
	V = T / 647.25e0;
	W = fabs(1.e0 - V);
	B = 0.e0;
	for (I = 1; I <= 8; I++)
	{
		L_Z = I;
		B = B + A[I] * pow(W, ((L_Z + 1.e0) / 2.e0));
	}
	Q = B / V;
	return (22.093e0 * exp(Q));
}

/* ---------------------------------------------------------------------- */
 LDBLE Phreeqc::
VLEST(LDBLE T)
/* ---------------------------------------------------------------------- */
/*
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	LDBLE A = -1.59259e1, B = 6.57886e-2, C = -1.12666e-4, D = 7.33191e-8,
		E = 1.60229e3, F = 2.88572e0, G = 650.0e0;

	return (A + B * T + C * T * T + D * T * T * T + E / T + F / (G - T));
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
DFIND(LDBLE * DOUT, LDBLE P, LDBLE D, LDBLE T)
/* ---------------------------------------------------------------------- */
/*
C
C     ROUTINE TO FIND DENSITY CORRESPONDING TO INPUT PRESSURE P(MPA), AN
C     TEMPERATURE T(K), USING INITIAL GUESS DENSITY D(G/CM3). THE OUTPUT
C     DENSITY IS IN G/CM3.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	int L;
	LDBLE DD, RT, PP_dfind, DPD, DPDX, DP, X;
	/*      LDBLE DD, RT, PP, DPD, DPDX, DP, X; */

	DD = D;
	RT = GASCON * T;
	if (DD <= 0.e0)
		DD = 1.e-8;
	if (DD > 1.9e0)
		DD = 1.9e0;
	L = 0;
	for (L = 1; L <= 30; L++)
	{
		if (DD <= 0.e0)
			DD = 1.e-8;
		if (DD > 1.9e0)
			DD = 1.9e0;
		QQ(T, DD);
		PP_dfind = RT * DD * BASE(DD) + Q0;
		DPD = RT * (Z + Y * DZ) + Q5;
		/*
		   C
		   C  THE FOLLOWING 3 LINES CHECK FOR NEGATIVE DP/DRHO, AND IF SO ASSUME
		   C    GUESS TO BE IN 2-PHASE REGION, AND CORRECT GUESS ACCORDINGLY.
		   C
		 */
		if (DPD <= 0.e0)
		{
			if (D >= .2967e0)
				DD = DD * 1.02e0;
			if (D < .2967e0)
				DD = DD * .98e0;
			if (L <= 10)
				continue;
		}
		else
		{
/* 13 */
			DPDX = DPD * 1.1e0;
			if (DPDX < 0.1e0)
				DPDX = 0.1e0;
			DP = fabs(1.e0 - PP_dfind / P);
			if (DP < 1.e-8)
				break;
			if (D > .3e0 && DP < 1.e-7)
				break;
			if (D > .7e0 && DP < 1.e-6)
				break;
			X = (P - PP_dfind) / DPDX;
			if (fabs(X) > .1e0)
				X = X * .1e0 / fabs(X);
			DD = DD + X;
			if (DD < 0.e0)
				DD = 1.e-8;
		}
	}
	if (L > 30)
		error_msg("In subroutine DFIND", STOP);
	//cite Jonathan Toner, remove error_msg
	//DD = 1;
/*   20 CONTINUE */
	*DOUT = DD;
	return OK;
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
QQ(LDBLE T, LDBLE D)
/* ---------------------------------------------------------------------- */
/*
C
C     THIS ROUTINE CALCULATES, FOR A GIVEN T(K) AND D(G/CM3), THE RESIDUL
C     CONTRIBUTIONS TO: PRESSURE (Q), DP/DRHO (Q5)
C     THIS SUBROUTINE IS USED IN DENSITY OF WATER CALCULATION.
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	/* COMMON /NCONST/ */
	LDBLE G[41] =
		{ 0, -.53062968529023e3, .22744901424408e4, .78779333020687e3,
		-.69830527374994e2, .17863832875422e5, -.39514731563338e5,
		.33803884280753e5,
		-.13855050202703e5, -.25637436613260e6, .48212575981415e6,
		-.34183016969660e6,
		.12223156417448e6, .11797433655832e7, -.21734810110373e7,
		.10829952168620e7,
		-.25441998064049e6, -.31377774947767e7, .52911910757704e7,
		-.13802577177877e7,
		-.25109914369001e6, .46561826115608e7, -.72752773275387e7,
		.41774246148294e6,
		.14016358244614e7, -.31555231392127e7, .47929666384584e7,
		.40912664781209e6,
		-.13626369388386e7, .69625220862664e6, -.10834900096447e7,
		-.22722827401688e6,
		.38365486000660e6, .68833257944332e4, .21757245522644e5,
		-.26627944829770e4,
		-.70730418082074e5, -.225e0, -1.68e0, .055e0, -93.0e0
	};
	/*int II[41]={0,4*0,4*1,4*2,4*3,4*4,4*5,4*6,4*8,2*2,0,4,3*2,4}; */
	int II[41] =
		{ 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5,
		5,
		5, 6, 6, 6, 6, 8, 8, 8, 8, 2, 2, 0, 4, 2, 2, 2, 4
	};
	/*int JJ[41]={0,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,1,3*4,0,2,0,0}; */
	int JJ[41] =
		{ 0, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3,
		5,
		7, 2, 3, 5, 7, 2, 3, 5, 7, 1, 4, 4, 4, 0, 2, 0, 0
	};
	int NC = 36;
	/* COMMON /ADDCON/ */
	LDBLE ATZ[5] = { 0, 64.e1, 64.e1, 641.6e0, 27.e1 }, ADZ[5] =
	{
	0, .319e0, .319e0, .319e0, 1.55e0}, AAT[5] =
	{
	0, 2.e4, 2.e4, 4.e4, 25.e0}, AAD[5] =
	{
	0, 34.e0, 4.e1, 3.e1, 1.05e3};
	LDBLE *QZT;
	LDBLE QR[12], QT[11] /*, QZT[10] */ ;
	/*EQUIVALENCE (QT(2),QZT(1)) */

	LDBLE E, Q10, Q20, V, QP, DDZ, DEL, EX1, DEX, ATT, TX,
		TAU, EX2, TEX, QM, FCT, Q5T;
	int I, K, L, J, KM;
	QZT = &(QT[1]);
	QR[1] = 0.e0;
	Q5 = 0.e0;
	Q0 = 0.e0;
	E = exp(-AA * D);
	Q10 = D * D * E;
	Q20 = 1.e0 - E;
	QR[2] = Q10;
	V = TZ / T;
	QT[1] = T / TZ;
	/*DO 4 I=2,10 */
	for (I = 2; I <= 10; I++)
	{
		QR[I + 1] = QR[I] * Q20;
		/*   4  QT[I]=QT[I-1]*V */
		QT[I] = QT[I - 1] * V;
	}
	/* DO 10 I=1,NC */
	for (I = 1; I <= NC; I++)
	{
		K = II[I] + 1;
		L = JJ[I];
		QP = G[I] * AA * QR[K + 1] * QZT[L];
		Q0 = Q0 + QP;
		/*10 Q5 = Q5 + AA*(2.e0/D-AA*(1.e0-E*(K-1)/Q20))*QP */
		Q5 = Q5 + AA * (2.e0 / D - AA * (1.e0 - E * (K - 1) / Q20)) * QP;
	}
	QP = 0.e0;
	/* DO 20 J=37,40 */
	for (J = 37; J <= 40; J++)
	{
		if (fabs(G[J]) < 1.0e-20)
			continue;
		K = II[J];
		KM = JJ[J];
		DDZ = ADZ[J - 36];
		DEL = D / DDZ - 1.e0;
		if (fabs(DEL) < 1.e-10)
			DEL = 1.e-10;
		EX1 = -AAD[J - 36] * pow(DEL, K);
		if (EX1 <= -88.028e0)
		{
			DEX = 0.e0;
		}
		else
		{
			DEX = exp(EX1) * pow(DEL, KM);
		}
		ATT = AAT[J - 36];
		TX = ATZ[J - 36];
		TAU = T / TX - 1.e0;
		EX2 = -ATT * TAU * TAU;
		if (EX2 <= -88.028e0)
		{
			TEX = 0.e0;
		}
		else
		{
			TEX = exp(EX2);
		}
		Q10 = DEX * TEX;
		QM = KM / DEL - K * AAD[J - 36] * pow(DEL, (K - 1));
		FCT = QM * D * D * Q10 / DDZ;
		Q5T =
			FCT * (2.e0 / D + QM / DDZ) - pow((D / DDZ),
											  2) * Q10 * (KM / DEL / DEL +
														  K * (K -
															   1) * AAD[J -
																		36] *
														  pow(DEL, (K - 2)));
		Q5 = Q5 + Q5T * G[J];
		QP = QP + G[J] * FCT;
		/*  20  CONTINUE */
	}
	Q0 = Q0 + QP;
	return OK;
}

/* ---------------------------------------------------------------------- */
 LDBLE Phreeqc::
BASE(LDBLE D)
/* ---------------------------------------------------------------------- */
/*
C
C     THIS FUNCTION CALCULATES THE Z (=PBASE/(DRT)) NEEDED FOR FUNCTION
C
C     FROM L. HAAR, J. S. GALLAGHER, AND G. S. KELL, (1984)
C
*/
{
	LDBLE X, Z0, DZ0;
/*
C
C     G1,G2 AND GF ARE THE ALPHA, BETA AND GAMMA FOR DENSITY OF WATER
C     CALCULATIONS.  B1 AND B2 ARE THE 'EXCLUDED VOLUME' AND '2ND VIRIAL
C     SUPPLIED BY THE SUBROUTINE BB(T), WHICH ALSO SUPPLIES THE 1ST AND
C     2ND DERIVATIVES WITH RESPECT TO T (B1T,B2T,B1TT,B2TT).
C
*/
	Y = .25e0 * B1 * D;
	X = 1.e0 - Y;
	Z0 = (1.e0 + G1 * Y + G2 * Y * Y) / pow(X, 3);
	Z = Z0 + 4.e0 * Y * (B2 / B1 - GF);
	DZ0 =
		(G1 + 2.e0 * G2 * Y) / pow(X,
								   3) + 3.e0 * (1.e0 + G1 * Y +
												G2 * Y * Y) / pow(X, 4);
	DZ = DZ0 + 4.e0 * (B2 / B1 - GF);
	return (Z);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
DC(LDBLE T)
/* ---------------------------------------------------------------------- */
/*
C
C     THIS FUNCTION CALCULATES THE RELATIVE DIELECTRIC CONSTANT AS A
C     FUNCTION OF TEMPERATURE, ASSUMING ONE ATMOSPHERE PRESSURE
C     ACCORDING TO D. J. BRADLEY AND K. S. PITZER, (1979)
C
*/
{
	LDBLE D1000, C, B;
	LDBLE U[10] = { 0, 3.4279e2, -5.0866e-3, 9.4690e-7, -2.0525e0, 3.1159e3,
		-1.8289e2, -8.0325e3, 4.2142e6, 2.1417e0
	};
	D1000 = U[1] * exp(U[2] * T + U[3] * T * T);
	C = U[4] + U[5] / (U[6] + T);
	B = U[7] + U[8] / T + U[9] * T;
	return (D1000 + C * log((B + VP * 1.01325e0) / (B + 1000.0e0)));
}
