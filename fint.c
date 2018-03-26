#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
/*#include <gsl/gsl_interp2d.h>*/
/*#include <gsl/gsl_spline2d.h>*/
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


double f (double s, void * params) {
  double alpha = *(double *) params;
  double f = log(alpha*s) / sqrt(s);
  return f;
}

int
main (void)
{
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  gsl_matrix* M = gsl_matrix_alloc(76,4);
  FILE* myfile = fopen("newtabletransfer.dat","r");
  gsl_matrix_fscanf(myfile, M);
  fclose(myfile);

  gsl_vector *v_T = gsl_vector_alloc (76);
  gsl_vector *v_n = gsl_vector_alloc (76);
  gsl_vector *v_k = gsl_vector_alloc (76);

  /* Carefull with this size */
  gsl_matrix_get_col(v_T, M, 3);
  gsl_matrix_get_col(v_n, M, 2);
  gsl_matrix_get_col(v_k, M, 0);

  double x_array[11], y_array[11], z_array[11];

  for (int j=0;j<11;j++) {
	x_array[j] = gsl_vector_get(v_n,j);
	y_array[j] = gsl_vector_get(v_k,j);
	z_array[j] = gsl_vector_get(v_T,j);	
	}
  
/*  gsl_interp_accel *kacc = gsl_interp_accel_alloc ();*/
  gsl_interp_accel *nacc = gsl_interp_accel_alloc ();
 /* gsl_spline2d *spline = gsl_spline2d_alloc (gsl_interp2d_bicubic, 196,196);

  gsl_spline2d_init(spline, x, y, z, 196, 196); 

 */
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline,10);

  gsl_spline_init(spline, x_array, z_array, 10);
/*
  for (xi = x[0]; xi < x[195]; xi += 1) {
	for (yi = y[0]; yi < y[195]; yi += 1) {	
        	zi = gsl_spline2d_eval(spline, xi, yi, nacc, kacc);
		printf ("%g %g %g \n", xi, yi, zi);
		}
	}
*/
  double x_test = 5.8e+16;
  double zi = gsl_spline_eval(spline, x_test, nacc);
  printf ("%f %.18f \n", x_test, zi);

  

  double result, error;
  double expected = -4.0;
  double alpha = 1.0;
  double omega_m = 0.308;
  double omega_r = 9.134e-5;
  double hub = 2.179e-18;
  double c = 3.0e8;
 

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
  
  for (int a=0;a<76;a++) {
	for(int b=0;b<4;b++) {
		printf ("%f ",gsl_matrix_get(M,a,b));
	}
	printf ("\n");
  }

  for (int p=0;p<76;p++) {
	printf ("%.18f ", gsl_vector_get(v_T,p));
	printf ("%f ", gsl_vector_get(v_n,p));
	printf ("%f ", gsl_vector_get(v_k,p));
	printf ("\n");
	}


  gsl_integration_workspace_free (w);
  gsl_matrix_free(M);
  gsl_vector_free(v_T);
  gsl_vector_free(v_n);
  gsl_vector_free(v_k); 
  gsl_spline_free(spline);
/*  gsl_spline2d_free(spline);*/
 /* gsl_interp_accel_free(kacc);*/
  gsl_interp_accel_free(nacc);

  return 0;
}
