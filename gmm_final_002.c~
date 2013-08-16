#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define REDUCTOR_PREFACTOR 6.0

#define  S 1026 /* 1025  ******** 256000 jeszcze pojdzie, potem core dump */ 

#define A_R (k_4*c_H) /* tak lub na odwrot, tym steruje przechodzac miedzy a_r i b_r */
#define B_R rho_0

#define STALA_DO_A_R 0.0 /* te dwa z programow do momentow, byl konflikt oznaczen */
#define STALA_DO_B_R 1.0

/* obecny program jest przerobka z:

automaton_007_bis_gsl.c

znajdujacego sie w katalogu GSL_ODE */

/*** historia ***/
/* 
Stara wersja automaton_002.c, ktora wczoraj (18 07 2012) sciagnalem z laptopa. Przerabiam ostroznie!

Usuwam tablice i dodaje wypisywanie "co ktorys"

Sprawdzam .....006.... dnia 02 08 2012, z notatkami analit. usprawniam (scalam sumy (petle for)).

z tego powstal automaton_007_bis_gsl.c
	       	
*/
/* WHAT IS NEW AS COMPARED TO automaton_007_bis_gsl.c ? */
/* co robie? przerabiam na wzor ssm_final_001.c, czyli dodaje wypisywanie parametrow modelu, d_0 do parmetrow */


/* PROBLEM!!!! */

/* this program is dedicated to the special case of the kernels, and not to the most general one. Needs to be generalized (later! not now!) */

/*****************************************************************************
*****************************************************************************
*****************************************************************************/

/*const int S = 4; */ /*alternatywa dla dyrektywy preprocesora powyzej */

const int N = 8000;  /* ilosc odcinkow dlugosci h dzielacych zadany odcinek b - a */

const int co_ile_krokow_wypisuje = 1; /* co ile  odcinkow dlugosci h nastepuje wypisanie wyniku */


double t1 = 800.0;  /*final time*/

int n;
int i, k, l, p, q; /* for the 'for' loops */

 
/**** auxilliary variables *********/

long double av_a1, av_M0_1, av_M1_1, av_M1_2, av_M2_1, av_M2_2;
long double av_xi_1_1, av_xi_1_2, av_xi_k_1, av_xi_k_2, av_xi_k_3, av_xi_k_4; 
long double av_xi_s_minus_1_1, av_xi_s_minus_1_2, av_xi_s_minus_1_3;

/* functions */


/******** < ... > **********/

long double mean(long double m_0, long double m_1, long double m_2);  
/* function returning mean cluster size */

long double std_deviation(long double m_0, long double m_1, long double m_2);  
/* function returning standard deviation of the cluster size */

/* the kernels */

long double K(int i_1, int i_2, long double kapp_0, long double kapp_1, long double kapp_2);

long double F(int i_1, int i_2);

long double R(int i_1, double c_HH, double kk_4, double rho_00);

long double A_K(int i_1, long double kapp_0, long double kapp_1, long double kapp_2);

long double B_K(int i_1, long double kapp_0, long double kapp_1, long double kapp_2);


/*parameters of the model. I add d_0 as in ssm_final_001.c */


struct rparams

	{
	const double c1;  /* kappa_0 */
	const double c2;  /* kappa_1 */
	const double c3;  /* kappa_2 */
	
	const double c4;  /* c_0 */
	const double c5;  /* c_H (concentration of  hydrazine, i.e., in general, the reducing agent */	

	const double c6;  /* k_1 */
	const double c7;  /* k_3 */
	const double c8;  /* k_4 */

	const double c9;  /* rho_0 */
	
	const double c10;  /* d_0 */ /*initial c_a, needed for a two-step WF */
	};


 



int func (double t, const double y[], double f[], void *params)

{
        double kappa_0 = (( struct rparams *) params) ->c1; 
	double kappa_1 = (( struct rparams *) params) ->c2; 
	double kappa_2 = (( struct rparams *) params) ->c3; 

	double c_0 	= (( struct rparams *) params) ->c4; 
	double c_H 	= (( struct rparams *) params) ->c5; 


	double k_1 = (( struct rparams *) params) ->c6; 
	double k_3 = (( struct rparams *) params) ->c7; 
	double k_4 = (( struct rparams *) params) ->c8; 

	double rho_0 = (( struct rparams *) params) ->c9;

	/* Dodane dopiero w tym programie, stad nie po kolei: */

	double d_0 	= (( struct rparams *) params) ->c10; /*initial c_a, needed for a two-step WF */ 

/* duze S tu to male s w rachunkach */
/* S + 4 equations for y[0] = c_p, y[i] = \xi_i, i = 1, ......., S-1, y[S] = M_0, y[S+1] = M_0, y[S+2] = M_2,  , y[S+3] = c_a */

/*********   c pi  or y[0] ******************************/


 f[0] =  - (k_1 * c_H) * y[0];


/********* y[1] or  xi_1  *******************************/

	av_xi_1_1 = 0.0;

		for(l = 1; l < S; l++)

		{
		av_xi_1_1 = av_xi_1_1 + K(1, l, kappa_0, kappa_1, kappa_2) * y[l];	
		}

	av_xi_1_2 = 0.0;

		for(p = 2; p < S; p++)

		{
		av_xi_1_2 = av_xi_1_2 + F(1, p-1) * y[p];	
		}




	f[1] =
			   ( (k_3 * c_H)  -  R(1, c_H,  k_4, rho_0) * y[1] ) * y[S+3]
				-  y[1] * (    av_xi_1_1 + ( A_K(1, kappa_0, kappa_1, kappa_2) * y[S+1] + B_K(1, kappa_0, kappa_1, kappa_2) * y[S] )     )
 				+  av_xi_1_2 ;

	

/*********  xi_k, k from  to s-2. !!! k = s-1 treated separately *******************************/


	av_xi_k_1 = 0.0; /* related to 0.5 * K[i][k-i] * y[i] * y[k-i]*/
	av_xi_k_2 = 0.0; /* related to - 0.5 * F[i][k-i]*/
	av_xi_k_3 = 0.0; /* related to -  K[k][l] * y[l] * y[k]   */
	av_xi_k_4 = 0.0;
 

 for(k = 2; k < S - 1; k++)

 {
 	av_xi_k_1 = 0.0; /* related to 0.5 * K[i][k-i] * y[i] * y[k-i]*/
	av_xi_k_2 = 0.0; /* related to - 0.5 * F[i][k-i]*/
	av_xi_k_3 = 0.0; /* related to -  K[k][l] * y[l] * y[k]   */
	av_xi_k_4 = 0.0;
	
	for(i = 1; i < k; i++)

		{
		av_xi_k_1 = av_xi_k_1 + K(i, k-i, kappa_0, kappa_1, kappa_2) * y[i] * y[k-i];	
	 	av_xi_k_2 = av_xi_k_2 + F(i, k-i);	
		}

	for(l = 1; l < S; l++)

		{
		av_xi_k_3 = av_xi_k_3 + K(l, k, kappa_0, kappa_1, kappa_2) * y[l];	
		}


		for(l = k+1; l < S; l++)

		{
		av_xi_k_4 = av_xi_k_4 + F(l-k, k) * y[l];	
		}


	f[k]   =   ( R(k-1, c_H,  k_4, rho_0) * y[k-1] - R(k, c_H,  k_4, rho_0) * y[k]) *y[S+3]
				+ 0.5 * av_xi_k_1 - 0.5 * y[k] * av_xi_k_2  
				-  y[k] * ( A_K(k, kappa_0, kappa_1, kappa_2) * y[S+1] + B_K(k, kappa_0, kappa_1, kappa_2) * y[S] + av_xi_k_3 ) 
				+ av_xi_k_4
			   ;

	 
 } /* for po k */

/*********  xi_k,  k = s - 1 *******************************/

	av_xi_s_minus_1_1 = 0.0; /* related to 0.5 * K[i][s_minus_1-i] * y[i] * y[s_minus_1-i]*/
	av_xi_s_minus_1_2 = 0.0; /* related to - 0.5 * F[i][s_minus_1-i]*/
	av_xi_s_minus_1_3 = 0.0; /* related to -  K[s_minus_1][l] * y[l] * y[s_minus_1]   */
	
	
	for(i = 1; i < S-1; i++)

		{
		av_xi_s_minus_1_1 = av_xi_s_minus_1_1 + K(i, S-1-i, kappa_0, kappa_1, kappa_2) * y[i] * y[S-1-i];	
		

		av_xi_s_minus_1_2 = av_xi_s_minus_1_2 + F(i, S-1-i);	
		}

	for(l = 1; l < S; l++)

		{
		av_xi_s_minus_1_3 = av_xi_s_minus_1_3 + K(l, S-1, kappa_0, kappa_1, kappa_2) * y[l];	
		}


	f[S-1]  =
			  ( R(S-2, c_H,  k_4, rho_0) * y[S-2] - R(S-1, c_H,  k_4, rho_0) * y[S-1]) *y[S+3]
			+ 0.5 * av_xi_s_minus_1_1 - 0.5 * av_xi_s_minus_1_2  * y[S-1]
			-  y[S-1] * ( A_K(S-1, kappa_0, kappa_1, kappa_2) * y[S+1] + B_K(S-1, kappa_0, kappa_1, kappa_2) * y[S] + av_xi_s_minus_1_3  ) 
			      ;


/*********  M_0 or  y[S] ******************************/

av_M0_1 = 0.0;

for(q = 1; q < S; q++)

	{

		for(i = q; i < S; i++)

		{
		av_M0_1 = av_M0_1 + K(i, S-1 + q-i, kappa_0, kappa_1, kappa_2) * y[i]* y[S-1 + q-i];	

		}

	}

f[S]  =  R(S-1, c_H,  k_4, rho_0) * y[S-1] *y[S+3] 
	+ 0.5 * av_M0_1 
	- 0.5 * kappa_0 * y[S] * y[S] - kappa_1 * y[S] * y[S+1] - 0.5 * kappa_2 * y[S+1] * y[S+1];
	

/********* M_1 or   y[S+1] ******************** pow(S-1.0 + q, 2.0) ***********/


av_M1_1 = 0.0;

for(q = 1; q < S; q++)

	{

		for(i = q; i < S; i++)
		{
		
	    av_M1_1 = av_M1_1 + (S-1.0 + (long double)q) * K(i, S-1 + q-i, kappa_0, kappa_1, kappa_2) * y[i]* y[S-1 + q-i];	

		}

	}

av_M1_2 = 0.0;

for(l = 1; l < S; l++)

		{
		av_M1_2 = av_M1_2 + (long double)l * y[l]* (A_K(l, kappa_0, kappa_1, kappa_2) * y[S+1] + B_K(l, kappa_0, kappa_1, kappa_2) * y[S]);	
		}

	

f[S+1] =
		 S * R(S-1, c_H,  k_4, rho_0) * y[S-1] * y[S+3] 
		 + y[S+3] * (A_R * y[S+1] + B_R * y[S]) 
		 + 0.5 * av_M1_1  
		 +  av_M1_2
		  ;



/********* M_2 or   y[S+2] ******************** pow(S-1.0 + q, 2.0) ***********/

av_M2_1 = 0.0;

for(q = 1; q < S; q++)

	{

		for(i = q; i < S; i++)
		{
		
	    av_M2_1 = av_M2_1 + pow(S-1.0 + q, 2.0) * K(i, S-1 + q-i, kappa_0, kappa_1, kappa_2) * y[i]* y[S-1 + q-i];	

		}

	}

av_M2_2 = 0.0;

for(l = 1; l < S; l++)

		{
		av_M2_2 = av_M2_2 + (long double)l * y[l] 
		* (
			2.0 * A_K(l, kappa_0, kappa_1, kappa_2) * y[S+2] 
			+ (2.0 * B_K(l, kappa_0, kappa_1, kappa_2) + (long double)l * A_K(l, kappa_0, kappa_1, kappa_2) ) * y[S+1] 
			+ (long double)l * B_K(l, kappa_0, kappa_1, kappa_2) * y[S]
     		  );	
		}




 f[S+2] =
		pow(S, 2.0) * R(S-1, c_H,  k_4, rho_0) * y[S-1] *y[S+3] 
		 +y[S+3] * (2.0 * A_R * y[S+2] + (2.0 * B_R + A_R) * y[S+1] + B_R * y[S]) 
		 + 0.5 * av_M2_1  
		 +  av_M2_2
		 +  kappa_0 * y[S+1] * y[S+1] + 2.0 * kappa_1 * y[S+1] * y[S+2] + kappa_2 * y[S+2] * y[S+2]
 		 ;




 
/*********   c alpha or y[S+3] ******************************/



av_a1 = 0.0;

	for(i = 1; i < S; i++)

	{
	av_a1 = av_a1 + R(i, c_H,  k_4, rho_0) * y[i];	

	}


 f[S+3] = (k_1 * c_H) * y[0] - (   (k_3 * c_H) + av_a1 + (A_R * y[S+1] + B_R * y[S])   ) * y[S+3];


 

return GSL_SUCCESS;
}






/***********************************************************************************************************
******************    M A I N     **********************************************************
******************************************************************************************/


int main(void) 
{

 
/***************** the begining of the GSL part **************************/

/*kappa_0, kappa_1, kappa_2,  c_0,                  c_H               k_1,     k_3,             k_4, 	                 rho_0 = \tilde{k}_4,                          d_0*/

struct rparams p = 

{0.0,   0.0,   0.0,   0.0*0.000075,    REDUCTOR_PREFACTOR * 0.001,   0.0*1.188,   0.0387,   STALA_DO_A_R *  7.53 * 100000,  STALA_DO_B_R * 7.53 * 100 * REDUCTOR_PREFACTOR, 0.000075};  

gsl_odeiv2_system sys = {func, NULL, S + 4, &p};     	   /*gsl_odeiv2_step_rkf45 */
gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-16, 1e-16, 0.0);

int i;
double t = 0.0;

double y[S + 4]; /* tablica zmiennych zaleznych (stezenia, momenty) */ 

long double stand_M_0;
long double stand_M_1;
long double stand_M_2;

long double aux_stand_M_0;
long double aux_stand_M_1;
long double aux_stand_M_2;



/* S + 4 variables: y[0] = c_p, y[i] = \xi_i, i = 1, ......., S-1, y[S] = M_0, y[S+1] = M_0, y[S+2] = M_2, y[S+3] = c_a */

 
/************* zerowanie reczne elem.  ************************/

for(i = 0; i < S + 4 ; i++)
	{
	y[i] = 0.0;
	 
/*	printf("# y[%d] = % .16Lf", i, y[i]); */
	}
 
    
/*******inicjalizacja tablicy zmiennych ******************/

y[0] = p.c4; /* c_p(0) = c_0 */
y[S+3] = p.c10; /* c_a(0) = d_0 */ 

 



printf ("# Values of the model parameters: \n#\n# Coagulation kernel parameters: kappa_0 = %.1f, kappa_1 = %.1f, kappa_2 = %.1f \n#\n# Precursors: c_0 = %.10f, d_0 = %.10f, Reductor: c_H = %.10f \n#\n# (exponential notation: c_0 = %.6e, d_0 = %.6e, c_H = %.6e)\n#\n# Reaction rate constants (bare!!): k_1 = %.10f, k_3 = %.10f, k_4 = %.10f, b_R = rho_0 = %.10f \n#\n# (exponential notation: k_1 = %.7e, k_3 = %.7e, k_4 = %.7e) \n#\n", p.c1, p.c2, p.c3, p.c4, p.c10, p.c5, p.c4, p.c10, p.c5, p.c6, p.c7, p.c8, p.c9, p.c6, p.c7, p.c8);  

printf ("# A_R = %.1f, B_R = %.1f \n#\n# N = %d, t1 = %.0f s \n#\n# co_ile_krokow_wypisuje = %d  \n#\n# s = %d\n#\n#\n#\n#\n#\n#\n", STALA_DO_A_R,  STALA_DO_B_R, N, t1, co_ile_krokow_wypisuje, S);

printf ("# State parameters: \n# \n#  t            c_p                c_a                 M_0                 M_1                  M_2               M_0^{S}	       M_1^{S}	          M_2^{S}        c_p + c_a + M_1      M_1/M_0     std_dev(size)           xi_1,             xi_2....\n#\n");


	for (i = 1; i <= N; i++) /*N = 100 in the orginal example */ 
	{

	double ti = i * t1 / N;
	int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

		if (status != GSL_SUCCESS)
		{
		printf ("error, return value=%d\n", status);
		break;
		}
	
		stand_M_0 = 0.0;
		stand_M_1 = 0.0;
		stand_M_2 = 0.0;

		aux_stand_M_0 = 0.0;
		aux_stand_M_1 = 0.0;
		aux_stand_M_2 = 0.0;

		for(k = 1; k < S ; k++)
		{
		aux_stand_M_0 = aux_stand_M_0 + y[k];	
		aux_stand_M_1 = aux_stand_M_1 + (long double)k * y[k];
		aux_stand_M_2 = aux_stand_M_2 + pow((long double)k, 2.0) * y[k];
	
		}
		stand_M_0 = aux_stand_M_0 + y[S];
		stand_M_1 = aux_stand_M_1 + y[S+1];
		stand_M_2 = aux_stand_M_2 + y[S+2];


		if (i % co_ile_krokow_wypisuje == 0)
		{

			printf("% .4f % .16f % .16f ", t, y[0], y[S+3]);
			printf("% .16Lf % .16Lf % .16Lf ", stand_M_0, stand_M_1, stand_M_2);
			printf("% .16f % .16f % .16f ", y[S], y[S+1], y[S+2]);
			printf("% .12Lf ", y[0] + y[S+3] + stand_M_1);	
			printf("% .12Lf ", mean(stand_M_0, stand_M_1, stand_M_2));
			printf("% .12Lf ", std_deviation(stand_M_0, stand_M_1, stand_M_2));
		/*		
			for(k = 1; k < S ; k++)
			{
				printf("% .16f ", y[k]);
			}
		*/

			for(k = 1; k < S ; k = 2 * k)  /*printf(" x[%d]=% .16f ", k, y[k]);*/
			{
				printf("% .16f ", y[k]);
			}
			printf(" \n"); 
	
	 

		} /* od 'if' od krokow czasowych */

		
		if (i == N)
		{
		printf("# \n");
		printf("# Teraz wypisuje asymptotyczne wartosci wszystkich stezen klastrow \n"); 
		for(k = 1; k < S ; k++)  /*printf(" x[%d]=% .16f ", k, y[k]);*/
			{
				printf("#7# %d  % .16f \n", k, y[k]);
			}
			printf("#M_0=#  % .16f \n", y[S]);
		}



	} /* for po 'i', glowna petla (kroki czasowe) */


	gsl_odeiv2_driver_free (d);


/***************** the end of the GSL part **************************/


printf ("# State parameters: \n# \n#  t            c_p                c_a                 M_0                 M_1                  M_2               M_0^{S}	       M_1^{S}	          M_2^{S}        c_p + c_a + M_1      M_1/M_0     std_dev(size)           xi_1,             xi_2....\n#\n");
 
  

 

return 0;
}  /*koniec funkcji main*/


/**** mean and standard deviation ***/


long double mean(long double m_0, long double m_1, long double m_2)
	{

	if (m_1 == 0.0 || m_0 == 0.0)
	return 0;
	else 
	return m_1/m_0;

	}

long double std_deviation(long double m_0, long double m_1, long double m_2)	{
	
	if (m_1 == 0.0 || m_0 == 0.0)
	return 0;
	else 	
	return  sqrt( m_2/m_0 - pow(m_1/m_0, 2.0));

	}

/****** Kernel functions ******************************/ 


long double K(int i_1, int i_2, long double kapp_0, long double kapp_1, long double kapp_2)
	{
	return kapp_2 * (long double)i_1 * (long double)i_2  + ( (long double)i_1 + (long double)i_2 ) * kapp_1  +   kapp_0;
	}


long double F(int i_1, int i_2)
	{
	return  0;
	}


long double R(int i_1, double c_HH, double kk_4, double rho_00)
	{
	return  c_HH * kk_4 * i_1 + rho_00;
	}


long double A_K(int i_1, long double kapp_0, long double kapp_1, long double kapp_2)
	{
	return i_1 * kapp_2  +  kapp_1;
	}


long double B_K(int i_1, long double kapp_0, long double kapp_1, long double kapp_2)
	{
	return i_1 * kapp_1  +  kapp_0;
	}
