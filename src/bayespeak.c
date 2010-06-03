#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include "bayespeak_gsl.h"
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>

#define debug_flag 0

//Global variables to aid debugging:
int globalstart, globalend;

FILE *my_file6;

//Look up functions

double ran0(const ran2_state_t* vstate)
{
	return ran2_get_double (vstate);
}

double normal(const ran2_state_t* vstate)
{
	return gsl_ran_gaussian_ratio_method (vstate, 1);
}

float normald(float x)
{
	//RIGHT TAIL
	return (1 - gsl_cdf_ugaussian_Q (x));
}

float factln(int n)
{
	return logfact(n);
}

//double lngammafn(double x) is OK

double gammarv(double ia, double ib, ran2_state_t* vstate)   
{
	return gsl_ran_gamma(vstate,ia,ib);
}

double betarv(double ia, double ib, ran2_state_t* vstate)  
{
	return gsl_ran_beta (vstate, ia, ib);
}

//Exit with error message

void errormsg(char msg[])
{
	error("Error in region %d-%d:\n%s", globalstart, globalend, msg);
}

void warnmsg(char msg[])
{
	warning("Warning in region %d-%d:\n%s", globalstart, globalend, msg);
}

//Print debug message

void say(char thing[])
{
	#if debug_flag
		Rprintf("%s",thing);
		fflush(stdout);
	#endif
}


// ***********************************
//		  Main algorithm
// ***********************************

//register routine with R
/*R_CallMethodDef callMethods[] = {*/
/*	{"bayespeak", &bayespeak, 13},*/
/*	{NULL, NULL, 0}*/
/*};*/

/*void R_init_myLib(DllInfo *info)*/
/*{*/
/*	R_registerRoutines(info, NULL, callMethods, NULL, NULL);*/
/*}*/


//int bayespeak(int *score2pos, int *score2neg, int *w, double *w_dbl, int *w_max, int *passtoN, char **chr, int *start, int *end, int *win, int *iterations, double *para_out)
SEXP bayespeak(SEXP score2posR, SEXP score2negR, SEXP wR, SEXP w_dblR, SEXP w_maxR, SEXP passtoNR, SEXP chrR, SEXP startR, SEXP endR, SEXP winR, SEXP iterationsR, SEXP para_outR, SEXP priorR) //FIXME
{
	//Conversions from R to C
	int* score2pos = INTEGER_POINTER(score2posR);
	int* score2neg = INTEGER_POINTER(score2negR);
	int* w = INTEGER_POINTER(wR);
	double* w_dbl = NUMERIC_POINTER(w_dblR);
	int* w_max = INTEGER_POINTER(w_maxR);
	int* passtoN = INTEGER_POINTER(passtoNR);
	//chrR doesn't need to be converted.
	int* start = INTEGER_POINTER(startR);
	int* end = INTEGER_POINTER(endR);
	int* win = INTEGER_POINTER(winR);
	int* iterations = INTEGER_POINTER(iterationsR);
	double* para_out = NUMERIC_POINTER(para_outR);
	double* prior = NUMERIC_POINTER(priorR);
	//end of conversions from R to C

	//globalise start/end for informative error
	int N = *passtoN; //= (*end-*start)/win;
	globalstart = *start;
	globalend = *end;

	int *stateZ;
	double *scorepos, *scoreneg, *probX;
	stateZ = (int*) R_alloc(N, sizeof(int));
	scorepos = (double*) R_alloc(N, sizeof(double));
	scoreneg = (double*) R_alloc(N, sizeof(double));
	probX = (double*) R_alloc(N, sizeof(double));

  ran2_state_t number;								//seed for GSL random number generator
  ran2_set(&number, 0);

  double rlamda[4];
  double sum1, sum2, sum3, sum4, sum5, sum5new, sum6, sum6new;
  sum1 = 0.0;
  sum3 = 0.0;
  double temp, temps[4];

  int iters = *iterations;
  double iters2 = (double)*iterations;
  int burnin = ceil(iters/2);
  int i, j, n, z, t, k1, k2;
  double  M1, M2, M3, M4, M, a, b;
  double u0, u1, u2, u3, u4, gamma, gammanew;
  double c0, c1, c2, cc0, cc1, cc2, z5, z6, z7;
  double aa0, aa1, C0, C1, C2, logA0, logA1, logA2;
  double Ap, Atheta, Bp, Btheta, Agamma, Bgamma;   // priors: p ~ Beta(Ap,Bp) and theta ~ Beta(Atheta,Btheta)
  double Aa0, Ba0, Ab0, Bb0, Aa1, Ba1, Ab1, Bb1;   // priors: a1 ~ Gamma(Aa1,Ba1) etc
  double loglhood;
  Ap = 2.0;
  Bp = 2.0;
  Atheta = 2.0;
  Btheta = 2.0;
  Agamma = 0.5;
  Bgamma = 1.0;
  //specifiable prior...
  Aa0 = prior[0]; //default: 5
  Ba0 = prior[1]; //default: 5
  Ab0 = prior[2]; //default: 10
  Bb0 = prior[3]; //default: 5
  Aa1 = prior[4]; //default: 25
  Ba1 = prior[5]; //default: 4
  Ab1 = prior[6]; //default: 0.5
  Bb1 = prior[7]; //default: 5

  double sd1, sd0, sd2;
  sd1 = 6.0;
  sd0 = 6.0;
  sd2 = 0.2;
  double theta, p, pi[4], rate[3];
  p = 0.20;
  theta = 0.90;
  pi[0] = 0.25;
  pi[1] = 0.25;
  pi[2] = 0.25;
  pi[3] = 0.25;
  rate[0] = 0.0;
  rate[1] = 0.0;
  rate[2] = 0.0;
  double lambda1, lambda0, a1, b1, a0, b0;
  a1=0.5;				 // starting points (can be changed) 
  b1=0.5;
  a0=0.5;
  b0=0.5;
  lambda0 = 5.0;
  lambda1 = 15.0;
  gamma = 1.0;

	int abort = 0;

	//initialise gamma_pow - lookup table to replace costly pow(gamma, w[i])
	double *gamma_pow = (double*) R_alloc((*w_max + 2), sizeof(double));
	temp = gamma;
	gamma_pow[0] = 1.0;
	gamma_pow[1] = gamma;
	for(i=2; i<= *w_max; i++)
	{
		temp *= gamma;
		gamma_pow[i] = temp;
	}
	
	double *gammanew_pow = (double*) R_alloc((*w_max + 2), sizeof(double));
	gammanew_pow[0] = 1.0;
	//pow(gamma, w[i]) will be equivalent to gamma_pow[w[i]].
	 
  double **Q, **Y, **cdfZ, **probZ, **P2;
  Q = (double**) R_alloc(4, sizeof(double*));
  Y = (double**) R_alloc(4, sizeof(double*));
  cdfZ = (double**) R_alloc(4, sizeof(double*));
  probZ = (double**) R_alloc(N, sizeof(double*));
  P2 = (double**) R_alloc(4, sizeof(double*));

	for(i=0; i<4; i++)
	{
		stateZ[i] = 0;
		Q[i] = (double*) R_alloc(4, sizeof(double));
  		for (j=0; j<4; j++)
		{
		   Q[i][j]=0.0;
		}
		Y[i] = (double*) R_alloc(N, sizeof(double));
		cdfZ[i] = (double*) R_alloc(N, sizeof(double));
		P2[i] = (double*) R_alloc(N, sizeof(double));
	}   
	for(i=0; i<N; i++)
	{ 
		probZ[i] = (double*) R_alloc(4, sizeof(double));
	}

/* reading in files in C would go here*/


	for (i=0; i<N; i++)
	{
		scorepos[i] = (double)score2pos[i];
		scoreneg[i] = (double)score2neg[i];
	}

	for (t=0; t<N-1; t++)
	{
		probX[t]=0.0;

		sum1 = sum1 + scorepos[t] + scoreneg[t+1];
		sum3 = sum3 + (scorepos[t] + scoreneg[t+1])*w_dbl[t];
	}
	probX[N-1] = 0.0;


	if(sum1 > FLT_MIN) //Go ahead only if there are reads. If not, forget it and output empty file.
	{
		for (n=0; n<(iters/10); n++)
		{
			for (z=0; z<10; z++)
			{
				u0 = ran0(&number);  
				u1 = ran0(&number);  
				u2 = ran0(&number);  
				u3 = ran0(&number);  
				M = 0.0;
				M1 = 0.0;
				M2 = 0.0;
				M3 = 0.0;
				M4 = 0.0;
				sum2 = 0.0;
				sum4 = 0.0;
				sum5 = 0.0;
				sum6 = 0.0;
				sum5new = 0.0;
				sum6new = 0.0;
				Q[0][0] = 1-p;
				Q[0][1] = p;
				Q[1][2] = 1-theta;
				Q[1][3] = theta;
				Q[2][0] = 1-p;
				Q[2][1] = p;
				Q[3][2] = 1-theta;
				Q[3][3] = theta;

				//calculate lhood of HMM
				//----------------------

				for (t=0; t<N-1; t++)
				{
					temp = gamma_pow[w[t]];
					sum5 = sum5 + temp;
					rlamda[0] = lambda0 * temp;
					rlamda[1] = (lambda0 + lambda1) * temp;
					rlamda[2] = rlamda[1];
					rlamda[3] = rlamda[1];

					k1 = score2pos[t];
					k2 = score2neg[t+1];
					for (i=0; i<4; i++)
					{
						Y[i][t] = exp((k1 + k2) * log(rlamda[i]) - rlamda[i] - rlamda[i] - factln(k1) - factln(k2));   // So the Y's are the probabilities of every set of counts for the underlying state Z_t = i
						//Y[i][t] = P(Y_t = data | Z_t = i)
					}
				}

				say("A");

				//get first column of P2, as probabilities proportional to vector pi * P(Z_0 = i)
				for (i=0; i<4; i++)
				{
					temps[i] = pi[i] * Y[i][0];
				}

				temp = 1.0/(temps[0]+temps[1]+temps[2]+temps[3]);
				for (i=0; i<4; i++)
				{
				  P2[i][0] = temps[i] * temp;
				}

				loglhood = log(P2[stateZ[0]][0]);

				say("B");

				//get subsequent columns of P2 (forward variable).
				for (t=1; t<N-1; t++)
				{
					for (i=0; i<4; i++)
					{
						temps[i] = (Q[0][i]*P2[0][t-1] + Q[1][i]*P2[1][t-1] + Q[2][i]*P2[2][t-1] + Q[3][i]*P2[3][t-1]) * Y[i][t];
					}

					temp = (temps[0]+temps[1]+temps[2]+temps[3]);
					if(temp < FLT_MIN) //if P2s will sum to 0...
					{
						//underflow
						//Rprintf("\nWarning: Underflow in P2 array, pass %d, at t = %d", 10*n + z + 1, t);
						for (i=0; i<4; i++)
						{
							P2[i][t] = 0.25;
						}				
					}
					else
					{
						//no underflow - update P2
						temp = 1.0/temp;

						for (i=0; i<4; i++)
						{
							P2[i][t] = temps[i] * temp;
						}
					}
				}
				cdfZ[0][N-2] = P2[0][N-2];  //cdfZ is CDF of Z
				probZ[N-2][0] = P2[0][N-2];

				if (u3<cdfZ[0][N-2])
				{
					stateZ[N-2] = 0;
				}
				else
				{
					for (i=1; i<4; i++)
					{
						//"cdfZ" only calculated if we need it
						cdfZ[i][N-2] = cdfZ[i-1][N-2] + P2[i][N-2];
						if (u3>cdfZ[i-1][N-2] && u3<cdfZ[i][N-2])
						{
						  stateZ[N-2] = i;
						}
					}
				}

				say("C");

				probX[N-2] = probZ[N-2][2] + probZ[N-2][3] + probX[N-2]; // post prob of X[t]=1
				probX[N-1] = probZ[N-2][1] + probZ[N-2][3] + probX[N-1];

				if (stateZ[N-2]>0)
				{
					sum2 = sum2 + scorepos[N-2] + scoreneg[N-1];
					sum4 = sum4 + (scorepos[N-2] + scoreneg[N-1]) * w_dbl[N-2];
					sum6 = sum6 + gamma_pow[w[N-2]];
				}
				  
				for (t=N-3; t>(-1); t--)
				{
					for (i=0; i<4; i++)
					{
						probZ[t][i] = P2[i][t] * Q[i][stateZ[t+1]];
					}

					temp = 1.0 / (probZ[t][0]+probZ[t][1]+probZ[t][2]+probZ[t][3]);

					for (i=0; i<4; i++)
					{
						probZ[t][i] = probZ[t][i] * temp;
					}

					cdfZ[0][t] = probZ[t][0];

					u4 = ran0(&number);  

					if (u4<cdfZ[0][t])
					{
						stateZ[t] = 0;
					}
					else
					{
						for (i=1; i<4; i++)
						{
							//"cdfZ" only calculated if we need it
							cdfZ[i][t] = cdfZ[i-1][t] + probZ[t][i];
							if (u4>cdfZ[i-1][t] && u4<cdfZ[i][t])
							{
								stateZ[t] = i;
							}
						}
					}
				}


				say("D");

				for (t=N-3; t>0; t--)
				{
					switch(stateZ[t])
					{
						case 0:
							M2++;
						break;
						case 1:
							M1++;
							sum2 += scorepos[t] + scoreneg[t+1];
							sum4 += (scorepos[t] + scoreneg[t+1]) * w_dbl[t];
							sum6 += gamma_pow[w[t]];
						break;
						case 2:
							M4++;
							sum2 += scorepos[t] + scoreneg[t+1];
							sum4 += (scorepos[t] + scoreneg[t+1]) * w_dbl[t];
							sum6 += gamma_pow[w[t]];
						break;
						case 3:
							M3++;
							sum2 += scorepos[t] + scoreneg[t+1];
							sum4 += (scorepos[t] + scoreneg[t+1]) * w_dbl[t];
							sum6 += gamma_pow[w[t]];
						break;
					}
				}


				say("E");
				  
				for (t=0; t<N-1; t++)
				{
					probX[t] = probZ[t][2] + probZ[t][3] + probX[t]; // post prob of X[t]=1
				}

				//end calculation of HMM lhood
				//-----------------------------

				say("F");

				//GIBBS SAMPLING update of p, theta.

				a = Ap + M1;	
				b = Bp + M2;

				do {
				p = betarv(a,b,&number);
				} while (p<=0.0001);

				a = Atheta + M4;
				b = Btheta + M3;

				do {
				theta = betarv(a,b,&number);		   
				} while (theta<=0.0001);

				//METROPOLIS-HASTINGS update of alpha_0, alpha_1, gamma
				//Proposal density is truncated normal distn. Truncation is to prevent alpha < 0
				//Probability of rejection is therefore at most half.

				do {
				z5 = normal(&number);
				aa0 = a0 + z5*sd0;									
				} while (aa0<=0.0001);

				do {
				z6 = normal(&number);
				aa1 = a1 + z6*sd1;		   
				} while (aa1<=0.0001);

				do {
				z7 = normal(&number);
				gammanew = gamma + z7*sd2;		   
				} while (gammanew<=0.0001);

				c0 = normald(-a0/sd0); 
				c1 = normald(-a1/sd1);
				c2 = normald(-gamma/sd2);
				cc0 = normald(-aa0/sd0); 
				cc1 = normald(-aa1/sd1);
				cc2 = normald(-gammanew/sd2);
				
				C0 = log((1.0-c0)/(1.0-cc0));
				C1 = log((1.0-c1)/(1.0-cc1));
				C2 = log((1.0-c2)/(1.0-cc2));

				logA1 = (Aa1 - 1.0) * log(aa1/a1) + (aa1-a1) * (log(lambda1) - Ba1 + log(b1)) + lngammafn(a1) - lngammafn(aa1) + C1;
				logA0 = (Aa0 - 1.0) * log(aa0/a0) + (aa0-a0) * (log(lambda0) - Ba0 + log(b0)) + lngammafn(a0) - lngammafn(aa0) + C0;

				say("G");

				//calculate gammanew_pow for use below:
				temp = gammanew;
				gammanew_pow[0] = 1.0;
				gammanew_pow[1] = gammanew;
				for(i=2; i<= *w_max; i++)
				{
					temp *= gammanew;
					gammanew_pow[i] = temp;
				}

				for (i=0; i<N-1; i++)
				{
					sum5new = sum5new + gammanew_pow[w[i]];
					if (stateZ[i]>0)
					{
						sum6new = sum6new + gammanew_pow[w[i]];
					}
				}
				say("H");
				logA2 = (sum3 + sum4 + Agamma - 1.0) * log(gammanew/gamma) + (gamma-gammanew) * Bgamma - 2 * lambda0 * (sum5new - sum5) - 2 * lambda1 * (sum6new - sum6) + C2;


				if (log(u0)<logA0)
				{
					 a0 = aa0;
					rate[0] = rate[0] + 1.0;
				}
				if (log(u1)<logA1)
				{
					 a1 = aa1;
					rate[1] = rate[1] + 1.0;
				}
				if (log(u2)<logA2)
				{
					gamma = gammanew;
					rate[2] = rate[2] + 1.0;

					//MH has set gamma = gammanew, so update gamma_pow = gammanew_pow
					for(i=*w_max; i >= 0; i--)
					{
						gamma_pow[i] = gammanew_pow[i];
					}
				}

				//End of METROPOLIS-HASTINGS
				//Gibbs update of b0, b1, then lambda0, lambda1 (already did a0, a1 with MH)

				int loopcount = 0;
				a = Ab0 + a0;
				b = 1.0/(Bb0 + lambda0);
				do {
					b0 = gammarv( a, b, &number); 
					loopcount++;
					if(loopcount > 100) {warnmsg("b0 too small - job aborted. "); abort = 1; break;}
				} while (b0<=0.001);
				if(abort){break;}

				loopcount = 0;
				a = Ab1 + a1;
				b = 1.0/(Bb1 + lambda1);
				do {
					b1 = gammarv( a, b, &number);
					loopcount++;
					if(loopcount > 100) {warnmsg("b1 too small - job aborted. "); abort = 1; break;}
				} while (b1<=0.001);
				if(abort){break;}

				//reset
				sum5 = 0.0;
				sum6 = 0.0;


				say("I");
				for (i=0; i<N-1; i++)
				{
					temp = gamma_pow[w[i]];
					sum5 += temp;

					if (stateZ[i]>0)
					{
						sum6 += temp;
					}
				}
				say("J");

				loopcount = 0;
				a = a0 + sum1;
				b = 1.0/(b0 + 2*sum5);
				do {
					lambda0 = gammarv( a, b, &number);
					loopcount++;
					if(loopcount > 100) {warnmsg("lambda0 too small - job aborted. "); abort = 1; break;}
				} while (lambda0<0.001);
				if(abort){break;}

				say("K");

				loopcount = 0;
				a = a1 + sum2;
				b = 1.0/(b1 + 2*sum6);
				do {
					lambda1 = gammarv( a, b, &number);   
					loopcount++;
					if(loopcount > 100) {warnmsg("lambda1 too small - job aborted. "); abort = 1; break;}
				} while (lambda1<0.001);
				if(abort){break;}

				say("L");

				//store model parameters for use in QC - reduce to mean later FIXME
				if(10*n + z > burnin)
				{
					para_out[0] += p;
					para_out[1] += theta;
					para_out[2] += lambda0;
					para_out[3] += lambda1;
					para_out[4] += gamma;
					para_out[5] += loglhood;
				}

			}
			R_CheckUserInterrupt(); //respond to Ctrl-C termination
		}

		Rprintf(".");
		fflush(stdout);
	}

	//finish taking mean of parameters (or kill them off if we aborted earlier)
	if(abort)
	{
		for(i=5; i>=0; i--)
		{
			para_out[i] = 0;
		}
	}
	else 
	{
		for(i=5; i>=0; i--)
		{
			para_out[i] = para_out[i]/(iters - burnin);
		}
	} //FIXME iters %% 10 != 0?

	//how much memory do we need?
	j = 0;
	if(!abort)
	{
		for (i=0; i<N-2; i++)
		{
			if (probX[i]>(0.01*iters2))
			{
				j++;
			}
		}
	}

	//collect enriched regions into R vectors
	SEXP outstart, outend, outPP;

	PROTECT(outstart = NEW_INTEGER(j));
	PROTECT(outend = NEW_INTEGER(j));
	PROTECT(outPP = NEW_NUMERIC(j));

	j = 0;
	if(!abort)
	{
		for (i=0; i<N-2; i++)
		{
			if (probX[i]>(0.01*iters2))
			{
				//chr?
				INTEGER_POINTER(outstart)[j] = *start+(i*(*win));
				INTEGER_POINTER(outend)[j] = *start+(i+1)*(*win);
				NUMERIC_POINTER(outPP)[j] = probX[i]/(iters2);
				j++;
			}
		}
	}

	//put regions and parameters into an output object "out"
	SEXP out, outnames;

	PROTECT(outnames = allocVector(STRSXP, 7));
	SET_STRING_ELT(outnames,0,mkChar("chr"));
	SET_STRING_ELT(outnames,1,mkChar("start"));
	SET_STRING_ELT(outnames,2,mkChar("end"));
	SET_STRING_ELT(outnames,3,mkChar("PP"));
	SET_STRING_ELT(outnames,4,mkChar("para"));
	SET_STRING_ELT(outnames,5,mkChar("jobstart"));
	SET_STRING_ELT(outnames,6,mkChar("jobend"));

	PROTECT(out = allocVector(VECSXP, 7));
	SET_VECTOR_ELT(out, 0, chrR);
	SET_VECTOR_ELT(out, 1, outstart);
	SET_VECTOR_ELT(out, 2, outend);
	SET_VECTOR_ELT(out, 3, outPP);
	SET_VECTOR_ELT(out, 4, para_outR);
	SET_VECTOR_ELT(out, 5, startR);
	SET_VECTOR_ELT(out, 6, endR);

	setAttrib(out, R_NamesSymbol, outnames);

	UNPROTECT(5);
	//end R objects


  return(out);
}

