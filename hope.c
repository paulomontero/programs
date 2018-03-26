#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>


/* Huge  run */
/*#define K_MAX 0.3
#define KPERPENDICULAR 0.005
#define NINT 600
#define DELTAK 0.0005
*/

/* maximum k perp run 
#define K_MAX 0.3
#define NINT 600
#define DELTAK 0.0005
#define KPERPENDICULAR 0.1
*/

/* What is going to be a large k perpendicular value? Maybe 0.1? ok. */

/* new big run, real upper plane
#define K_MAX 0.03
#define KPERPENDICULAR 0.005
#define NINT 60
#define DELTAK 0.0005
*/

/*NEW BIG RUN WITH KPERP biggest = 0.1, upper plane*2 scheme*/
#define K_MAX 0.3
#define KPERPENDICULAR 0.0010
#define NINT 600
#define DELTAK 0.0005






int
main (void)
{

/* GRAB, INTERPOLATE AND INTERPOLATE 2D THE CLASS OUPUT FROM TRANSFER.PY */
  gsl_matrix* M1 = gsl_matrix_alloc(7570635,4);
  printf ("Matrix worked \n");
  FILE* myfile = fopen("newtabletransfer.dat","r");
  gsl_matrix_fscanf(myfile, M1);
  printf ("Matrix was scanned \n");
  fclose(myfile);
  
  long int sum;
  sum =0;
  long int i;
  long int a;
  a = 0;
  static double t_array[1939*2000];
  double n_array[1939];
  for (long int j=0;j<7570634;j++) {
  	if (gsl_matrix_get(M1,j,0) < gsl_matrix_get(M1,j+1,0)) {
		i = i + 1;
		double xi_array[i], zi_array[i];
		for (long int b=0;b<i;b++) {	
			xi_array[b] = gsl_matrix_get(M1,b+sum,2);
			zi_array[b] = gsl_matrix_get(M1,b+sum,3);
 			}
		gsl_interp_accel *xacc = gsl_interp_accel_alloc ();
		gsl_spline *splinei = gsl_spline_alloc (gsl_interp_cspline,i);
  		gsl_spline_init(splinei, xi_array, zi_array, i);
		for(int s=0; s<1939;s++){  					
			n_array[s] = 5.83e+16 + s*5.1468e+13;
			t_array[s+1939*a] = gsl_spline_eval(splinei,n_array[s],xacc);
			}
		sum = sum + i;
		a = a + 1;
		i = 0;
  		gsl_interp_accel_free(xacc);
  		gsl_spline_free(splinei);
	} else {
		i = i + 1; 
  	}
	}
  gsl_matrix_free(M1);
  double k_array[2000];

  for(int p=0; p<2000;p++){
  	k_array[p] = 0.0001 + p*0.0005;
	}

  gsl_interp_accel *kacc = gsl_interp_accel_alloc ();
  gsl_interp_accel *nacc = gsl_interp_accel_alloc ();
  gsl_spline2d *spline = gsl_spline2d_alloc (gsl_interp2d_bicubic, 1939, 1999);

  gsl_spline2d_init(spline, n_array, k_array, t_array, 1939, 1999); 


  
  /* Here I start the actual integration part */

/* IMPORTANT PARAMETERS */
  /*double omega_m = 0.308000000;
  double omega_r = 9.134000000e-5;
  double hub = 2.17900000e-18;*/
  double c = (299792.458000) / (3.086e+19); 


 
  printf ("Careful with units, speed of light is in Mpc/s = %0.18f \n", c);

  double etasurfacels = 5.83000000e+16;
/* TESTS FOR GSL INTEGRATION */

  

 

  gsl_matrix* M2 = gsl_matrix_alloc(1999*1999,5);
  printf ("Matrix worked \n");
  FILE* myfile2 = fopen("LAtabla2000.dat","r");
  gsl_matrix_fscanf(myfile2, M2);
  printf ("Matrix was scanned \n");
  fclose(myfile2);
  
  static double fabs1[1999*1999],freal[1999*1999], fimag[1999*1999];
  double kmag[1999], kpar[1999];
  for (int b2 = 0; b2 < 2000; b2++) {
	kmag[b2] = 0.0001 + b2*0.0005;
	kpar[b2] = 0.0001 + b2*0.0005;
  }


  for (int j2 = 0; j2 < 1999*1999; j2++) {	
	fabs1[j2] = gsl_matrix_get(M2,j2,2);
	freal[j2] = gsl_matrix_get(M2,j2,3);
	fimag[j2] = gsl_matrix_get(M2,j2,4);
  }

  gsl_matrix_free(M2);
  printf("uhmmm %e %e %e\n", freal[1999*1999-1], fimag[1999*1999-1], fabs1[1999*1999-1]);

/* HERE I INTERPOLATE THE FABS */
  gsl_interp_accel *kmagacc2 = gsl_interp_accel_alloc ();
  gsl_interp_accel *kparacc2 = gsl_interp_accel_alloc ();
  gsl_spline2d *splinefabs2 = gsl_spline2d_alloc (gsl_interp2d_bicubic, 1999, 1999);

  gsl_interp_accel *kmagacc3 = gsl_interp_accel_alloc ();
  gsl_interp_accel *kparacc3 = gsl_interp_accel_alloc ();
  gsl_spline2d *splinefre = gsl_spline2d_alloc (gsl_interp2d_bicubic, 1999, 1999);


  gsl_interp_accel *kmagacc4 = gsl_interp_accel_alloc ();
  gsl_interp_accel *kparacc4 = gsl_interp_accel_alloc ();
  gsl_spline2d *splinefimag = gsl_spline2d_alloc (gsl_interp2d_bicubic, 1999, 1999);
  
  gsl_spline2d_init(splinefabs2, kmag, kpar, fabs1, 1999, 1999);
  gsl_spline2d_init(splinefre, kmag, kpar, freal, 1999, 1999);
  gsl_spline2d_init(splinefimag, kmag, kpar, fimag, 1999, 1999);
  


/* Primordial power spectrum stuff */
/* Grab the A_s constant from transfer.py and n_s */
  double A_s = 2.206e-9;
  double n_s = 0.9652;
  double k_pivot = 0.05;

  double inline primordial_pk (double k_in) {
/* this is in units of Mpc^3 */
    return A_s*(pow(k_in/k_pivot,n_s-1))*(pow(k_in,-3));
  }

  	 
  printf ("Did it. \n");


/*should try to grab stuff that is common everywhere like the sin^4 in another way instead that in every single integrand, have to be careful with this. */
/* This is the common factor in all integrals from the spherical harmonics, it is everywhere */
  double sin4theta1sin4theta2 (double k1per, double k1m, double k2per, double k2m){
    return pow(k1per/k1m,4)*pow(k2per/k2m,4); 
  } 

/* HERE I COMPUTE stuff need for terms like j1, m1, k1, l1 */  
  double fabs_sq (double magk1, double k1par1) {
    return gsl_spline2d_eval(splinefabs2, magk1, k1par1, kmagacc2, kparacc2);
  }
/* Here is the factor fre(k1)fre(k2)-Fima(k1)Fima(k2) this comes in j2, k2 */ 
  double realffactor (double magk1, double k1par1, double magk2, double k2par1) {
    return (gsl_spline2d_eval(splinefre, magk1, k1par1, kmagacc2, kparacc2)*gsl_spline2d_eval(splinefre, magk2, k2par1, kmagacc2, kparacc2) - gsl_spline2d_eval(splinefimag, magk1, k1par1, kmagacc2, kparacc2)*gsl_spline2d_eval(splinefimag, magk2, k2par1, kmagacc2, kparacc2));
  }

/* Here is the factor fre(k1)fima(k2)+Fima(k1)Fre(k2) this is needed in l2, m2 */ 
  double imagffactor (double magk1, double k1par1, double magk2, double k2par1) {
    return (gsl_spline2d_eval(splinefre, magk1, k1par1, kmagacc2, kparacc2)*gsl_spline2d_eval(splinefimag, magk2, k2par1, kmagacc2, kparacc2) + gsl_spline2d_eval(splinefimag, magk1, k1par1, kmagacc2, kparacc2)*gsl_spline2d_eval(splinefre, magk2, k2par1, kmagacc2, kparacc2));
  }

/* cos (4 phi_1) */
  double cos4phi1 (double k1x, double k1per) {
    return 8*pow(k1x/k1per,4) - 8*pow(k1x/k1per,2) + 1;
  }

/* sin(4 phi_1) */
  double sin4phi1 (double k1x, double k1y, double k1per) {
    return (8*k1y*pow(k1x,3))/(pow(k1per,4)) - 4*(k1x*k1y)/(k1per*k1per);
  }

/* sin(4 phi_2) for the ones with positive k2 -> + 2 k1x kperp*/
  double sin4phi2pos (double k1x, double k1y, double k2perp) {
    return (8*k1y*pow(k1x+KPERPENDICULAR,3))/(pow(k2perp,4)) - (4*(k1x + KPERPENDICULAR)*k1y)/(k2perp*k2perp);
  }

/* cos (4 phi_2) for positive k2*/
  double cos4phi2pos (double k1x, double k2perp) {
    return 8*pow((k1x + KPERPENDICULAR)/(k2perp),4) - 8*pow((k1x+ KPERPENDICULAR)/(k2perp),2) + 1;
  }




/* Now let's define a kmax and make it change */

/* First kmax going to be 0.01 this implies N =20 */
/* Second kmax is 0.03 this implies N = 60 */

/* cfactor and dfactor comes from defining (cos4phi2 - i sin4phi2)(cos4phi1 + i sin4phi1)= C + i D,
   SSimilarly for (cos4phi1 - i sin4phi1)(cos4phi2 + i sin4phi2) = E + i F, note that E != than C since E has the negative values of 2k1xKEPERPENDICULAR. Seems like F andD vanish because integrand is odd. Instead of even. CHECK THIS ARGUMENT LATER */
  double result_j1, k1par1_a, k2par1_a, k1x_a, k1y_a, sin4sfactor_possign, result_m1, realffactor_possign, imagffactor_possign, Tek2_possign, Tek1, kmagni1, kmagni2_possign, result_k2re, result_k2im, result_l2re, result_l2im, cfactor, result_k1re, result_l1re, result_j2re, result_j2im, result_m2re, result_m2im; 
  double delta = pow(DELTAK,4);
  double k1per, k2per_possign, fabsfactork1, fabsfactork2;


 /* for (int ii1y = -1*NINT; ii1y < NINT; ii1y++) { this is for all*/
/*for upper plane we have */
  for (int ii1y = 0; ii1y < NINT; ii1y++) {
/* for lower plane
  for (int ii1y = -1*NINT; ii1y < 0; ii1y++) {*/
/* This for runing over upper	k1y_a = 0.00025 + DELTAK*ii1y;*/
/*	k1y_a = -0.00975 + DELTAK*ii1y;*/
/*	k1y_a = -0.29975 + DELTAK*ii1y;*/
	k1y_a = (ii1y + 0.5)*DELTAK;
	printf("Vamos por y loop: %d \n", ii1y);

  	for (int ii1x = -1*NINT; ii1x < NINT; ii1x++) {/*This for NINT 20 */
	/*	k1x_a = -0.00975 + DELTAK*ii1x;*/
	/*	k1x_a = -0.29975 + DELTAK*ii1x;*/
		k1x_a = (ii1x + 0.5)*DELTAK;
		/*This is for NINT 600 k1x_a = - 0.29975 + DELTAK*ii1x;*/
		/*printf("Vamos por x loop: %d \n", ii1x);*/
		/* PERPENDIICULAR KS*/
		k1per = sqrt(k1x_a*k1x_a + k1y_a*k1y_a);
		k2per_possign = sqrt(k1per*k1per + KPERPENDICULAR*KPERPENDICULAR + 2*k1x_a*KPERPENDICULAR);
		
		if (k1per < K_MAX && k2per_possign < K_MAX){
			for (int ii1p = 0; ii1p < NINT; ii1p++) {
	/*			k1par1_a = 0.0001 + DELTAK*ii1p; */
				k1par1_a = (ii1p + 0.5)*DELTAK;
				kmagni1 = sqrt(k1per*k1per + k1par1_a*k1par1_a);
				Tek1 = gsl_spline2d_eval(spline, etasurfacels, kmagni1, nacc, kacc);
			/* SPHERICAL HARMONICS FACTOR DEALING EXPLICTLY WITH PHI STUFF*/
				cfactor = cos4phi2pos(k1x_a, k2per_possign)*cos4phi1(k1x_a, k1per) + sin4phi2pos(k1x_a, k1y_a, k2per_possign)*sin4phi1(k1x_a, k1y_a, k1per);
			/* d and f terms vanish once integrated over ky */

				fabsfactork1 = fabs_sq(kmagni1,k1par1_a);
				for (int ii2p = 0; ii2p < NINT; ii2p++) {
		/*		k2par1_a = 0.0001 + DELTAK*ii2p;*/
					k2par1_a = (ii2p + 0.5)*DELTAK;
					kmagni2_possign = sqrt(k2par1_a*k2par1_a + k2per_possign*k2per_possign);
					if (kmagni1 < K_MAX && kmagni2_possign < K_MAX) {
				/*if (kmagni2_possign < K_MAX) {*/
				/*if (kmagni1 < K_MAX) { */
				/*if (k1per < K_MAX) { this one also works better, what's happening? */
				/*if (k2per_possign < K_MAX) {*/
					/*T = 1;*/
					/*printf("estoy en possigns %f %f %f %f \n", kmagni1, kmagni2_possign, lala, K_MAX);*/
						sin4sfactor_possign = sin4theta1sin4theta2(k1per, kmagni1, k2per_possign, kmagni2_possign)*primordial_pk(kmagni1)*primordial_pk(kmagni2_possign); 
							
						realffactor_possign = realffactor(kmagni1, k1par1_a, kmagni2_possign, k2par1_a);
						imagffactor_possign = imagffactor(kmagni1, k1par1_a, kmagni2_possign, k2par1_a);
						fabsfactork2 = fabs_sq(kmagni2_possign,k2par1_a);

						Tek2_possign = gsl_spline2d_eval(spline, etasurfacels, kmagni2_possign, nacc, kacc);
					
					/*if(kmagni1 > 0.9*K_MAX) {
						T = (K_MAX - kmagni1)/(0.1*K_MAX);
					}*/
						result_j1 += fabsfactork1*Tek2_possign*Tek2_possign*sin4sfactor_possign*(delta);						
				
						result_k2re += realffactor_possign*Tek1*Tek2_possign*sin4sfactor_possign*delta;
			
						result_k2im += imagffactor_possign*Tek1*Tek2_possign*sin4sfactor_possign*delta;
						result_k1re += fabsfactork1*Tek2_possign*Tek2_possign*sin4sfactor_possign*(cfactor)*delta;
						result_j2re += realffactor_possign*cfactor*Tek1*Tek2_possign*sin4sfactor_possign*(delta);/* - imagffactor_possign*dfactor*Tek1*Tek2_possign*sin4sfactor_possign*delta);*/

						result_j2im += imagffactor_possign*cfactor*Tek1*Tek2_possign*sin4sfactor_possign*delta; /* + realffactor_possign*dfactor*Tek1*Tek2_possign*sin4sfactor_possign*delta);*/

						result_m1 += fabsfactork2*Tek1*Tek1*sin4sfactor_possign*delta;
						result_l2re += realffactor_possign*Tek1*Tek2_possign*sin4sfactor_possign*delta;
						result_l2im += -1*imagffactor_possign*Tek1*Tek2_possign*sin4sfactor_possign*delta;
						result_l1re += fabsfactork2*Tek1*Tek1*cfactor*sin4sfactor_possign*delta;
						result_m2re += realffactor_possign*Tek1*Tek2_possign*cfactor*sin4sfactor_possign*delta;
						result_m2im += -1*imagffactor_possign*Tek1*Tek2_possign*cfactor*sin4sfactor_possign*delta;
					} else {
					/*printf("breaking %d \n", ii2p);*/
						break; 
					}
				}
			}
		}
	}
  }

		

  printf("the  kmax and kperp for upper*2 plane run: %f %f \n", K_MAX, KPERPENDICULAR);
  printf("La integral j1: %e \n", result_j1*2);
  
  printf("La integral m1: %e \n", result_m1*2);
  printf("La integral k2 da : %e + i %e \n", result_k2re*2, result_k2im*2); 
  printf("La integral l2 da : %e + i %e \n", result_l2re*2, result_l2im*2);
  printf("La integral k1 da : %e + i 0 \n", result_k1re*2);
  printf("La integral l1 da : %e + i 0 \n", result_l1re*2); 
  printf("La integral j2 da : %e + i %e \n", result_j2re*2, result_j2im*2);   
  printf("La integral m2 da : %e + i %e \n", result_m2re*2, result_m2im*2);
 
  
  double finalresultreal, finalresultimag;
  finalresultreal =  result_j1*2 + result_m1*2 - result_k2re*2 - result_l2re*2 - result_k1re*2 - result_l1re*2 + result_j2re*2 +  result_m2re*2;
  finalresultimag = result_j2im*2 + result_m2im*2 - result_k2im*2 - result_l2im*2;
  printf("Entonces j1+j2-k1-k2-l1-l2+m1+m2 : %e + i %e \n", finalresultreal, finalresultimag);
 
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(kacc);
  gsl_interp_accel_free(nacc);
  gsl_interp_accel_free(kmagacc2);
  gsl_interp_accel_free(kparacc2);
  gsl_spline2d_free(splinefabs2);
  gsl_interp_accel_free(kmagacc3);
  gsl_interp_accel_free(kparacc3);
  gsl_spline2d_free(splinefre);
  gsl_interp_accel_free(kmagacc4);
  gsl_interp_accel_free(kparacc4);
  gsl_spline2d_free(splinefimag);
  return 0;
}
