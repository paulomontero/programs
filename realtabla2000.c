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

#define K_MAX 0.03

double f (double y, void * params) {
  double alpha = *(double *) params;
  double f = log(alpha*y) / sqrt(y);
  return f;
}

double etafactor (double yval, double om, double or, double h) {
  double factor = ((yval*yval)*(h*h)*(om*om) - 4*or) / (4*om);
  return factor;
} 

double complex expfactor (double yvar, double kvar, double cs, double etalss) {
  double xfact = cs*kvar*(yvar - etalss);
  double complex zvar = 0 + xfact * I;
  double complex expo = cexp(zvar); 
  return expo;
}


int
main (void)
{


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
/*		double *xi, *zi;
		int sizhow big can an index in a loop be in ce = i;	
		xi = (double*)calloc(size, sizeof(int));			
		zi = (double*)calloc(size, sizeof(int));	*/
		for (long int b=0;b<i;b++) {	
			xi_array[b] = gsl_matrix_get(M1,b+sum,2);
			zi_array[b] = gsl_matrix_get(M1,b+sum,3);
		/*	printf ("Esta: %f \n", xi_array[b]);*/
 			}
		printf ("Copiando datos \n");
		int n = sizeof(xi_array) / sizeof(xi_array[0]);
		printf ("Size of x array: %d \n", n);
		gsl_interp_accel *xacc = gsl_interp_accel_alloc ();
		gsl_spline *splinei = gsl_spline_alloc (gsl_interp_cspline,i);
  		gsl_spline_init(splinei, xi_array, zi_array, i);
		printf ("Interpolando \n");
		for(int s=0; s<1939;s++){  			
			/*printf ("El ciclo esta bien \n");*/			
			n_array[s] = 5.83e+16 + s*5.1468e+13;
			/*printf ("Logre hacer este n para, %d %f \n", s, n_array[s]);*/
			t_array[s+1939*a] = gsl_spline_eval(splinei,n_array[s],xacc);
			/*printf ("Logre hacer este t para %d \n", s);*/
			}
		sum = sum + i;
		a = a + 1;
		printf ("Sum: %ld \n", sum);
		printf ("Borrando todo \n");
		i = 0;
  		gsl_interp_accel_free(xacc);
  		gsl_spline_free(splinei);
	} else {
		i = i + 1; 
	/*	printf ("Contando filas \n");*/
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

  double n_test = 8.91e+16;
  double k_test = 0.0112;
  double t_test = gsl_spline2d_eval(spline, n_test, k_test, nacc, kacc);
  printf ("Esta: %f %f %0.18f \n", n_test, k_test, t_test);


  double n_test1 = 1e+17;
  double k_test1 = 0.0213;
  double t_test1 = gsl_spline2d_eval(spline, n_test1, k_test1, nacc, kacc);
  printf ("Esta: %f %f %0.18f \n", n_test1, k_test1, t_test1);
  
  /* Here I start the actual integration part */

  double omega_m = 0.308000000;
  double omega_r = 9.134000000e-5;
  double hub = 2.17900000e-18;
  double c = (299792.458000) / (3.086e+19); 

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;
  double expected = -4.0; 
  double alpha = 1.0;

  gsl_function F;
  F.function = &f;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        w, &result, &error); 

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);

  double etafac = etafactor (n_test1, omega_m, omega_r, hub);
  printf ("Factor de eta = %.18f \n", etafac);
 
  printf ("Careful with units, speed of light is in Mpc/s = %0.18f \n", c);

  double etasurfacels = 5.83000000e+16;
  double kpruebapara = 0.02;
  double kpruebaper = 0.009;
  double complex efac = expfactor (n_test1, kpruebapara, c,  etasurfacels);
  printf ("Exponencial, parte real, parte imaginaria = %f, %f \n", creal(efac), cimag(efac));

  struct my_f_params { double om; double or; double h; double nlss; double kgivenpara; double kgivenmag; double cs; };

  double my_f (double eta, void * p) {
    struct my_f_params * params = (struct my_f_params *)p;
    double om = (params->om);
    double or = (params->or);
    double h = (params->h);
    double nlss = (params->nlss);
    double kgivenpara = (params->kgivenpara);
    double cs = (params->cs);
    double kgivenmag = (params->kgivenmag);

    return cs*gsl_spline2d_eval(spline, eta, kgivenmag, nacc, kacc)*(pow(4*om / ((eta*eta)*(h*h)*(om*om) - 4 *or),4))*cos(cs*kgivenpara*(eta - nlss));
  }


  double my_f_imaginary (double eta, void * p) {
    struct my_f_params * params = (struct my_f_params *)p;
    double om = (params->om);
    double or = (params->or);
    double h = (params->h);
    double nlss = (params->nlss);
    double kgivenpara = (params->kgivenpara);
    double cs = (params->cs);
    double kgivenmag = (params->kgivenmag);

    return cs*gsl_spline2d_eval(spline, eta, kgivenmag, nacc, kacc)*(pow(4*om / ((eta*eta)*(h*h)*(om*om) - 4 *or),4))*sin(cs*kgivenpara*(eta - nlss));
  }

  gsl_function G_imaginary;
  gsl_function G;
  struct my_f_params params = { omega_m, omega_r, hub, etasurfacels, kpruebapara, kpruebaper, c};
  

  G.function = &my_f;
  G.params = &params;
  printf ("Integrand: %g \n", GSL_FN_EVAL (&G,n_test1));

  G_imaginary.function = &my_f_imaginary;
  G_imaginary.params = &params;
  printf ("Integrand: %g \n", GSL_FN_EVAL (&G_imaginary,n_test1));
  

  gsl_integration_workspace * u = gsl_integration_workspace_alloc (100000);


  gsl_integration_workspace * u_imaginary = gsl_integration_workspace_alloc (100000);

  double intresult, interror;
  double intresultima, interrorima;

  gsl_integration_qag (&G, etasurfacels, 15.73e+16, 1e-4, 1e-4, 100000, 6, u, &intresult, &interror);
  
  gsl_integration_qag (&G_imaginary, etasurfacels, 15.73e+16, 1e-4, 1e-4, 100000, 6, u_imaginary, &intresultima, &interrorima);


  /* This is working, but for ow it only have one of the integrals and only the real part!!!! */
  /*gsl_integration_qags (&G, etasurfacels, 15.73e+16, 0, 1e-7, 10000,
                        u, &intresult, &interror); */

  
  printf ("result  real part       = % .18f\n", intresult);
  printf ("estimated error = % .18f\n", interror);
  printf ("intervals       = %zu\n", u->size);


  
  printf ("result imaginary part         = % .18f\n", intresultima);
  printf ("estimated error = % .18f\n", interrorima);
  printf ("intervals       = %zu\n", u_imaginary->size);


  /* Now let's try to print a table of |F(k1//,ki_)|^2 */
/*  
  int ii = 0;
  int jj = 0;
  double k1paraj1_array[331], k1perpj1_array[331], f1real_array[331], f1imag_array[331], f2real_array[331], f2imag_array[331], fabs_array[331], ftotreal_array[331], ftotimag_array[331];
  int k0; */
/*  FILE *f1abs = fopen("f1abs","w");*/
/*
  for(int kp = 0; kp <= 9000; kp = kp + 500) {
	
	double kpar = kp*1.000000000000000000000/10000;
	if ( 1*1 - kp*kp < 0) {
		k0 = 0;
	} else {
		k0 = sqrt(1 - kp*kp); 
  		}
	for(int kpp = k0; kpp < sqrt(9400*9400 - kp*kp); kpp = kpp +500) {	
		double kper = kpp*1.000000000000000000000/10000;
		struct my_f_params params = { omega_m, omega_r, hub, etasurfacels, kpar, kper, c};
		printf ("Integrating for = % .18f %.18f \n", kper, kpar);		
		G.params = &params;
		G_imaginary.params = &params;
		gsl_integration_qag (&G, etasurfacels, 15.73e+16, 1e-3, 1e-3, 100000, 6, u, &intresult, &interror);
  		gsl_integration_qag (&G_imaginary, etasurfacels, 15.73e+16, 1e-3, 1e-3, 100000, 6, u_imaginary, &intresultima, &interrorima);
		k1paraj1_array[ii] = kpar;
		k1perpj1_array[ii] = kper;
		f1real_array[ii]   = intresult;
		f1imag_array[ii]   = intresultima;
		f2real_array[ii]   = intresult;
		f2imag_array[ii]   = intresultima; 
		fabs_array[ii]     = f1real_array[ii]*f1real_array[ii] + f1imag_array[ii]*f1imag_array[ii];
		ftotreal_array[ii] = f1real_array[ii]*f2real_array[ii] - f1imag_array[ii]*f2imag_array[ii];
		ftotimag_array[ii] = f1real_array[ii]*f2imag_array[ii] + f2real_array[ii]*f1imag_array[ii];
		printf(" %0.18f  %0.18f  %f %f %f\n", kpar, kper, fabs_array[ii], ftotreal_array[ii], ftotimag_array[ii]);*/
		/*fprintf(f1abs, " %0.18f  %0.18f  %f \n", kpar, kper, intresult*intresult + intresultima*intresultima); 		*/
			/*
		ii = ii + 1;		
		}
	ii = ii + 1;	
	}
*/
/* Here the new approach */
/*  size_t neval; */
  

  double kpartest = 0.2236;
  double kmagtest = 0.0331;
  size_t nputos = 430;
  gsl_integration_glfixed_table * aqui_imag = gsl_integration_glfixed_table_alloc(nputos);
  G_imaginary.params = &params;
  double result1 = gsl_integration_glfixed(&G_imaginary, etasurfacels, 15.7000e+16, aqui_imag);

  printf ("La integral para 430 puntos da: %f \n", result1);
  gsl_integration_glfixed_table_free(aqui_imag);

  nputos = 800;
  gsl_integration_glfixed_table * aqui_imag2 = gsl_integration_glfixed_table_alloc(nputos);
 
  double result2 = gsl_integration_glfixed(&G_imaginary, etasurfacels, 15.7000e+16, aqui_imag2);

  printf ("La integral para 800 puntos da: %f \n", result2);
  gsl_integration_glfixed_table_free(aqui_imag2);

  nputos = 30;
  gsl_integration_glfixed_table * aqui_imag3 = gsl_integration_glfixed_table_alloc(nputos);
 
  double result3 = gsl_integration_glfixed(&G_imaginary, etasurfacels, 15.7000e+16, aqui_imag3);

  printf ("La integral para 30 puntos da: %f \n", result3);
  gsl_integration_glfixed_table_free(aqui_imag3);  

  nputos = 1921;
  gsl_integration_glfixed_table * aqui_imag4 = gsl_integration_glfixed_table_alloc(nputos);
 
  double result4 = gsl_integration_glfixed(&G_imaginary, etasurfacels, 15.7000e+16, aqui_imag4);

  printf ("La integral para 1921 puntos da: %f \n", result4);
  gsl_integration_glfixed_table_free(aqui_imag4);
  
  printf ("Las diferencias con la de 1921 puntos dan: %f %f %f  \n", (result4-result1)*100/result4, (result4-result2)*100/result4, (result4-result3)*100/result4);

 
 /* FILE *tabla;
  tabla = fopen("LAtabla2000.dat","w");*/
/* Time to make the cycle and store the values on the arrays */
/* HERE EVERYTHING SHOULD BE 2000 */
 /* double the_kmag, the_kpar;
  double the_fabs, the_freal, the_fimag;

  double delta_eta, N_points, intresult_real, intresult_imag;
  for (int ii = 0; ii < 1999; ii++) {
	the_kmag = 0.000100000+ ii*0.00050000;*/
/* I am going to change stuff for testing, currently is a little bit slow for < 2000 both and change of 0.0005*ii same for jj to 0.05 which implies  < 20 */
/*	for (int jj = 0; jj < 1999; jj++) { 
		the_kpar = 0.00010000 + jj*0.00050000;
		if(the_kpar < 0.005) {
			delta_eta = 1/(2*the_kpar*c*100);
			N_points = (15.7e+16 - etasurfacels)/delta_eta + 1;
			nputos = round(N_points);
		} else {
			delta_eta = 1/(2*the_kpar*c);
			N_points = (15.7e+16 - etasurfacels)/delta_eta + 1;
			nputos = round(N_points);
		}			
		gsl_integration_glfixed_table * tabla_real = gsl_integration_glfixed_table_alloc(nputos);
		gsl_integration_glfixed_table * tabla_imag = gsl_integration_glfixed_table_alloc(nputos);
	
		struct my_f_params params = { omega_m, omega_r, hub, etasurfacels, the_kpar, the_kmag, c};
		printf ("Integrating for kmag, kpar = % .18f %.18f \n", the_kmag,  the_kpar);		
		G.params = &params;
		G_imaginary.params = &params;
		intresult_real = gsl_integration_glfixed(&G, etasurfacels, 15.7000e+16, tabla_real);
		intresult_imag = gsl_integration_glfixed(&G_imaginary, etasurfacels, 15.7000e+16, tabla_imag);

			
		the_fabs = intresult_real*intresult_real + intresult_imag*intresult_imag;
		the_freal = intresult_real;
		the_fimag = intresult_imag;
		fprintf(tabla, "%0.18f %0.18f %e %e %e \n", the_kmag, the_kpar, the_fabs, the_freal, the_fimag);

			
		
		gsl_integration_glfixed_table_free(tabla_real);
		gsl_integration_glfixed_table_free(tabla_imag);

	}
  }
  fclose(tabla);
*/

  gsl_matrix* M2 = gsl_matrix_alloc(1999*1999,5);
  printf ("Matrix worked \n");
  FILE* myfile2 = fopen("LAtabla2000.dat","r");
  gsl_matrix_fscanf(myfile2, M2);
  printf ("Matrix was scanned \n");
  fclose(myfile2);
  
  static double fabs3[1999*1999],freal3[1999*1999], fimag3[1999*1999];
  double kmag[1999], kpar[1999];
  for (int b2 = 0; b2 < 2000; b2++) {
	kmag[b2] = 0.0001 + b2*0.0005;
	kpar[b2] = 0.0001 + b2*0.0005;
  }

  for (int j2 = 0; j2 < 1999*1999; j2++) {	
	fabs3[j2] = gsl_matrix_get(M2,j2,2);
	freal3[j2] = gsl_matrix_get(M2,j2,3);
	fimag3[j2] = gsl_matrix_get(M2,j2,4);
  }


 /* printf("Last row fabs %.10e \n", fabs3[398]);
  printf("Last row fabs %.10e \n", fabs3[399]);
  printf("Last row MATRIX %.10e \n", gsl_matrix_get(M2,399,2));

  printf("Last row MATRIX %.10e \n", gsl_matrix_get(M2,398,2));

  
  printf("First row fabs %.10e \n", fabs3[0]);
  printf("First row MATRIX %.10e \n", gsl_matrix_get(M2,0,2));
  int da = sizeof(fabs3) / sizeof(fabs3[0]);
  printf ("Size of x array: %d \n", da);*/
  gsl_matrix_free(M2);
  gsl_interp_accel *kmagacc = gsl_interp_accel_alloc ();
  gsl_interp_accel *kparacc = gsl_interp_accel_alloc ();
  gsl_spline2d *splinefabs = gsl_spline2d_alloc (gsl_interp2d_bicubic, 1999, 1999);

  gsl_spline2d_init(splinefabs, kmag, kpar, fabs3, 1999, 1999);

  double kmagprueba1 = 0.250;
  double kparprueba1= 0.0003;
  double fabsprueba1 = gsl_spline2d_eval(splinefabs, kmagprueba1, kparprueba1, kmagacc, kparacc);
  printf ("Esta: %f %f %f \n", kmagprueba1, kparprueba1, fabsprueba1);



  
  double n_testa = 1e+17;
  double k_testa = 0.9991;
  double t_testa = gsl_spline2d_eval(spline, n_testa, k_testa, nacc, kacc);
  printf ("Esta: %f %f %e \n", n_testa, k_testa, t_testa);

  gsl_integration_workspace_free (w);
  gsl_integration_workspace_free (u);
  gsl_integration_workspace_free (u_imaginary);
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(kacc);
  gsl_interp_accel_free(nacc);

  gsl_interp_accel_free(kmagacc);
  gsl_interp_accel_free(kparacc);
  gsl_spline2d_free(splinefabs);
  return 0;
  }
/* Ran perfectly for the change of 20 instead of 2000 and change of 0.05 */
/* Also ran perfectly for 200 and change of 0.005, probably took like 30 mins */

/* Now I am goint to interpolat the stuff for fabs_matrix[ii][jj] since I need it */
