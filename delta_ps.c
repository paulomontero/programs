#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include <complex.h>

/*
  USAGE: delta_ps <deltax filename> <xH filename> <output filename> <redshift of box>

  box is assumed to be of the HII dimension defined in ANAL_PARAM.H

  P: This program will be modified to take two boxes instead of  one, with the objective of getting P_{m,xhi}(k)
*/

#define FORMAT (int) 0 /* 0= unpadded binary box; 1= FFT padded binary box (outdated) */
#define CONVERT_TO_DELTA (int) 0 /* 1= convert the field to a zero-mean delta; 
		     be careful not to do this with fields which are already zero mean */

int main(int argc, char ** argv){
  char filename[100];
  FILE *F;
  float REDSHIFT;
  
  int x,y,z, format;
  fftwf_complex *deltax;
  //P: NEED ANOTHER ONE OF THIS FFTWF
  fftwf_complex *deltaxH;
  //P: DO I NEED ANOTHER PLAN??? HOW DOES THE PLAN IN FFTWF WORKS?
  fftwf_plan plan;
  fftwf_plan planH;
  float k_x, k_y, k_z, k_mag, k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor;
  int i,j,k, n_x, n_y, n_z, NUM_BINS;
  double dvdx, ave, aveH, new_ave, new_aveH, *p_box, *pimag_box, *k_ave, *ptest_box;
  unsigned long long ct, *in_bin_ct;

  // check arguments
  // P: ADDING ANOTHER ARGUMENT HERE FOR THE NEW BOX.
  if (argc != 5){
    fprintf(stderr, "USAGE: delta_ps <deltax filename> <xH filename> <output filename> <redshift of box>\nAborting\n");
    return -1;
  }
  // initialize and allocate thread info
  if (fftwf_init_threads()==0){
    fprintf(stderr, "init: ERROR: problem initializing fftwf threads\nAborting\n.");
    return -1;
  }
  fftwf_plan_with_nthreads(NUMCORES); // use all processors for init

  //P: TWO AVERAGE TEMPS
  ave=0;
  aveH=0;

  //P: REDSHIFT OF BOX FOR CONVENIENCE
  double redz =  strtod(argv[4], NULL);

  //allocate and read-in the density array
  //P: I NEED ANOTHER ONE OF THIS.
  deltax = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  deltaxH = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  //P: I GUESS I CAN ALSO KEEP THE ERRORS FOR BOTH BOXES.

  if (!deltax){
    fprintf(stderr, "delta_T: Error allocating memory for deltax box\nAborting...\n");
    fftwf_cleanup_threads(); return -1;
  }
  if (!deltaxH){
    fprintf(stderr, "delta_T?: Error allocating memory for deltam box\nAborting...\n");
    fftwf_cleanup_threads(); return -1;
  }    
  F = fopen(argv[1], "rb");
  switch (FORMAT){
    // FFT format
  case 1:
    fprintf(stderr, "Reading in FFT padded deltax box\n");
    if (mod_fread(deltax, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS, 1, F)!=1){
      fftwf_free(deltax);
      fprintf(stderr, "deltax_ps.c: unable to read-in file\nAborting\n");
      fftwf_cleanup_threads(); return -1;
    }
    break;

    // unpaded format
  case 0:
    fprintf(stderr, "Reading in unpadded box\n");
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if (fread((float *)deltax + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	    fprintf(stderr, "init.c: Read error occured!\n");
	    fftwf_free(deltax);
	    fftwf_cleanup_threads(); return -1;	    
	  }
      	  ave += *((float *)deltax + HII_R_FFT_INDEX(i,j,k));
	}
      }
    }
    ave /= (double)HII_TOT_NUM_PIXELS;
    fprintf(stderr, "Average is %e\n", ave);
    break;

  default:
    fprintf(stderr, "Wrong format code\naborting...\n");
    fftwf_free(deltax);
    fftwf_cleanup_threads(); return -1;	    
  }
  fclose(F);

  //P: NEED SOMETHING SIMILAR TO THE PREVIOUS LOOP TO GET DELTAXH, HOWEVER, WHERE IN THE BOX IS fluctuation XH???

  F = fopen(argv[2], "rb");
  switch (FORMAT){
    // FFT format
  case 1:
    fprintf(stderr, "Reading in FFT padded deltax box\n");
    if (mod_fread(deltaxH, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS, 1, F)!=1){
      fftwf_free(deltaxH);
      fprintf(stderr, "deltax_ps.c: unable to read-in file\nAborting\n");
      fftwf_cleanup_threads(); return -1;
    }
    break;

    // unpaded format
  case 0:
    fprintf(stderr, "Reading in unpadded box\n");
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  if (fread((float *)deltaxH + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	    fprintf(stderr, "init.c: Read error occured!\n");
	    fftwf_free(deltaxH);
	    fftwf_cleanup_threads(); return -1;	    
	  }
	//P:NEED ANOTHER AVERAGE!
      	  aveH += *((float *)deltaxH + HII_R_FFT_INDEX(i,j,k));
	}
      }
    }
    aveH /= (double)HII_TOT_NUM_PIXELS;
    fprintf(stderr, "Average is %e\n", aveH);
    break;

  default:
    fprintf(stderr, "Wrong format code\naborting...\n");
    fftwf_free(deltaxH);
    fftwf_cleanup_threads(); return -1;	    
  }
  fclose(F);

  if (CONVERT_TO_DELTA){
    new_ave = 0;
    new_aveH = 0;
    fprintf(stderr, "Now converting field to zero-mean delta\n");
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  *((float *)deltax + HII_R_FFT_INDEX(i,j,k)) /= ave;
	  *((float *)deltax + HII_R_FFT_INDEX(i,j,k)) -= 1;
	  new_ave += *((float *)deltax + HII_R_FFT_INDEX(i,j,k));

	  *((float *)deltaxH + HII_R_FFT_INDEX(i,j,k)) /= aveH;
	  *((float *)deltaxH + HII_R_FFT_INDEX(i,j,k)) -= 1;
	  new_aveH += *((float *)deltaxH + HII_R_FFT_INDEX(i,j,k));
	}
      }
    }
    new_ave /= (double) HII_TOT_NUM_PIXELS;
    fprintf(stderr, "The mean value of the field is now %e\n", new_ave);
    new_aveH /= (double) HII_TOT_NUM_PIXELS;
    fprintf(stderr, "The mean value of the field is now %e\n", new_aveH);
  }


  fprintf(stderr, "Made it to before the FFTs\n");

  // do the FFTs
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax, (fftwf_complex *)deltax, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fprintf(stderr, "Did the plan for density\n");

  //P: PLAN FOR H FOURIER STUFF
  planH = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltaxH, (fftwf_complex *)deltaxH, FFTW_ESTIMATE);
  fftwf_execute(planH);
  fftwf_destroy_plan(planH);
  fprintf(stderr, "Did the plan for xh\n");
  fftwf_cleanup();

  fprintf(stderr, "Clean up and ready to loop through fftws\n");

  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
     deltax[ct] *= VOLUME/(HII_TOT_NUM_PIXELS+0.0);
     deltaxH[ct] *= VOLUME/(HII_TOT_NUM_PIXELS+0.0);
  }


  /******  PRINT OUT THE POWERSPECTRUM  *********/

  k_factor = 1.4;
  k_first_bin_ceil = DELTA_K;
  k_max = DELTA_K*HII_DIM;
  // initialize arrays
   NUM_BINS = 0;
  k_floor = 0;
  k_ceil = k_first_bin_ceil;
  while (k_ceil < k_max){
    NUM_BINS++;
    k_floor=k_ceil;
    k_ceil*=k_factor;
  }

  fprintf(stderr, "Initialized the arrays for power\n");

  p_box =  (double *)malloc(sizeof(double)*NUM_BINS);
  pimag_box = (double *)malloc(sizeof(double)*NUM_BINS);
  ptest_box = (double *)malloc(sizeof(double)*NUM_BINS);
  k_ave =  (double *)malloc(sizeof(double)*NUM_BINS);
  in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*NUM_BINS);

  if (!p_box || !in_bin_ct || !k_ave){ // a bit sloppy, but whatever..
    fprintf(stderr, "delta_T.c: Error allocating memory.\nAborting...\n");
    fftwf_free(deltax);
    fftwf_free(deltaxH);
    fftwf_cleanup_threads(); return -1;
  }

  fprintf(stderr, "Allocated memory for power\n");

  for (ct=0; ct<NUM_BINS; ct++){
    p_box[ct] = k_ave[ct] = 0;
    pimag_box[ct] = 0;
    in_bin_ct[ct] = 0;
  }


  // now construct the power spectrum file
  for (n_x=0; n_x<HII_DIM; n_x++){
    if (n_x>HII_MIDDLE)
      k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
    else
      k_x = n_x * DELTA_K;

    for (n_y=0; n_y<HII_DIM; n_y++){
      if (n_y>HII_MIDDLE)
	k_y =(n_y-HII_DIM) * DELTA_K;
      else
	k_y = n_y * DELTA_K;

      for (n_z=0; n_z<=HII_MIDDLE; n_z++){ 
	k_z = n_z * DELTA_K;
	
	k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

	// now go through the k bins and update
	ct = 0;
	k_floor = 0;
	k_ceil = k_first_bin_ceil;

	while (k_ceil < k_max){
	  // check if we fall in this bin
	  if ((k_mag>=k_floor) && (k_mag < k_ceil)){
	    in_bin_ct[ct]++;
 	    //P: NEED TO FIND OUT HOW TO MULTIPLY BY COMPLEX ONE
	    //P: STORE REAL PART IN PBOX, DO NOT STORE IMAGINARY PART BECAUSE INTEGRATING OVER HALF

	    p_box[ct] += pow(k_mag,3)*(creal(deltax[HII_C_INDEX(n_x, n_y, n_z)])*creal(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]) + cimag(deltax[HII_C_INDEX(n_x, n_y, n_z)])*cimag(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]))/(2.0*PI*PI*VOLUME); 
	    
	    pimag_box[ct] += pow(k_mag,3)*pow(cabs(deltax[HII_C_INDEX(n_x, n_y, n_z)]), 2) / (2.0*PI*PI*VOLUME);

	    //p_box[ct] += pow(k_mag,3)*(creal(deltax[HII_C_INDEX(n_x, n_y, n_z)])*creal(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]) + cimag(deltax[HII_C_INDEX(n_x, n_y, n_z)])*cimag(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]) + I*(cimag(deltax[HII_C_INDEX(n_x, n_y, n_z)])*creal(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]) - creal(deltax[HII_C_INDEX(n_x, n_y, n_z)])*cimag(deltaxH[HII_C_INDEX(n_x, n_y, n_z)])))/(2.0*PI*PI*VOLUME);
	    
	    //p_box[ct] += pow(k_mag,3)*(creal(deltax[HII_C_INDEX(n_x, n_y, n_z)])*creal(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]) + cimag(deltax[HII_C_INDEX(n_x, n_y, n_z)])*cimag(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]))/(2.0*PI*PI*VOLUME);


            //pimag_box[ct] += pow(k_mag,3)*(cimag(deltax[HII_C_INDEX(n_x, n_y, n_z)])*creal(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]) - creal(deltax[HII_C_INDEX(n_x, n_y, n_z)])*cimag(deltaxH[HII_C_INDEX(n_x, n_y, n_z)]))/(2.0*PI*PI*VOLUME);
	   
	    //p_box[ct] +=  pow(k_mag,3)*pow(cabs(deltax[HII_C_INDEX(n_x, n_y, n_z)]), 2) / (2.0*PI*PI*VOLUME); -> P: THIS IS THE ORIGINAL

	    //pimag_box[ct] += pow(k_mag,3)*(deltax[0]*deltax[0] + deltax[1]*deltax[1]) / (2.0*PI*PI*VOLUME); -> P: THIS DID NOT WORK

	    //ptest_box[ct] += pow(k_mag,3)*(creal(deltax[HII_C_INDEX(n_x,n_y,n_z)])*creal(deltax[HII_C_INDEX(n_x,n_y,n_z)]) + cimag(deltax[HII_C_INDEX(n_x,n_y,n_z)])*cimag(deltax[HII_C_INDEX(n_x,n_y,n_z)])) / (2.0*PI*PI*VOLUME); -> P: THE TEST WAS PERFECT

	    // note the 1/VOLUME factor, which turns this into a power density in k-space

	    //P: I THOUGHT BECAUSE OF PARITY THIS CROSS POWER SPECTRUM WAS SUPPOSED TO BE REAL
	    //P: QUITE WEIRD STUFF, MAKE SURE TO RUN WITH BOXES OF SAME REDSHIFT

            k_ave[ct] += k_mag;
	    break;
	  }

	  ct++;
	  k_floor=k_ceil;
	  k_ceil*=k_factor;

	}
      }
    }
  } // end looping through k box
  fftwf_free(deltax);
  fftwf_free(deltaxH);

  // now lets print out the k bins
  F = fopen(argv[3], "w");
  if (!F){
    fprintf(stderr, "delta_T.c: Couldn't open file %s for writting!\n", filename);
    fftwf_cleanup_threads(); return -1;
  }

  fprintf(stderr, "Preparing to print final po\n");
  for (ct=1; ct<NUM_BINS; ct++){
 //   fprintf(F, "%e\t%e\t%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), pimag_box[ct]/(in_bin_ct[ct]+0.0), ptest_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
    fprintf(F, "%e\t%e\t%e\t%e\t%e\n", redz, k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), pimag_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
  }
  fclose(F);

  /****** END POWER SPECTRUM STUFF   ************/

  free(p_box); free(k_ave); free(in_bin_ct); free(pimag_box); free(ptest_box);

  fftwf_cleanup_threads(); return 0;
}
