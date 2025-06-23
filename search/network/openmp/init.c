#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <getopt.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

#include "auxi.h"
#include "struct.h"
#include "init.h"
#include "settings.h"
#include "iniparser.h"

#if defined(_OPENMP)
#include <omp.h>
#endif



void read_ini_file( Search_settings *sett, 
		    Command_line_opts *opts,
		    int argc, 
		    char* argv[]) {
  
  char ini_fname[FNAME_LENGTH];
  dictionary *ini;
  int error = 0;
  
  if (argc > 1) {
    strcpy (ini_fname, argv[1]);
  } else if (argc > 2) {
    printf("WARNING: too many arguments, only first one is used.\n");    
  } else {
    printf("ERROR: missing input file name. Call: \"<executable> search.ini\" \n");
    exit(EXIT_FAILURE);
  }

  printf ("Loading config file %s\n", ini_fname);
  if ((ini = iniparser_load(ini_fname)) == NULL) {
    perror(ini_fname);
    exit(EXIT_FAILURE);
  }
  
  iniparser_dump(ini, stdout);
  
  // directory containing input data
  opts->indir = iniparser_getstring(ini, "search:indir", NULL);
  // output directory
  opts->outdir = iniparser_getstring(ini, "search:outdir", ".");
  // band number
  opts->band = iniparser_getint(ini, "search:band", 0);
  // time segment number
  opts->seg = iniparser_getint(ini, "search:seg", 0);
  // hemisphere 1,2 or 0 for both
  opts->hemi = iniparser_getint(ini, "search:hemi", 0);
  // candidate fstatistic threshold
  opts->thr = iniparser_getdouble(ini, "search:thr", 20.);
  // number of days (integer)
  sett->nod = iniparser_getint(ini, "search:nod", 0);
  // sampling interval of the input time series [seconds] (double)
  sett->dt = iniparser_getdouble(ini, "search:dt", -1.);
  // bands overlap [0-1.] defines band base frequency, fpo
  opts->overlap = iniparser_getdouble(ini, "search:overlap", 0.);
  // [0-0.5] range of frequencies to output; 0 is nothing, 0.5 is full band
  opts->narrowdown = iniparser_getdouble(ini, "search:narrowdown", -1.);
  // use data from subset of detectors only (default is to use all available)
  opts->usedet = iniparser_getstring(ini, "search:usedet", "");
  // name of the range file to read
  opts->range_file = iniparser_getstring(ini, "search:range_file", "");
  // name of the range file to dump and exit
  opts->dump_range_file = iniparser_getstring(ini, "search:dump_range_file", "");
  // name of the file with signall to add
  // TODO: single line signals to enable multiple signals, more flexible gsize
  opts->addsig = iniparser_getstring(ini, "search:addsig", "");
  // optional label of input and output files
  opts->label = iniparser_getstring(ini, "search:label", "");  
  
  // fstatistics normalization; if NULL white noise is assumed
  // currently only NULL and old blocks_avg are implemented
  opts->fstat_norm = iniparser_getstring(ini, "search:fstat_norm", NULL);

  // directed search options
  // right ascension 
  opts->ra = iniparser_getstring(ini, "search:ra", ""); // right ascension in radians [0, 2*PI]
  // declination
  opts->dec = iniparser_getstring(ini, "search:dec", ""); // declination in radians [-PI/2, PI/2]

  //   FLAGS
  // [0, 1] (in future other values can select different subsets of lines
  opts->veto_flag = iniparser_getint(ini, "search:veto_flag", 0);
  // [0, 1] if 1 generate vlines file and exit
  opts->gen_vlines_flag = iniparser_getint(ini, "search:gen_vlines_flag", 0);
  // [0, 1] generate checkpoint files at each new sky position
  opts->checkp_flag = iniparser_getint(ini, "search:checkp_flag", 0);

  
  // state_file is set in read_checkpoint

  
  // various checks
  if (! opts->indir) {
    error = 1; printf("[ERROR] missing indir !\n");
  }
  if (opts->band < 1) {
    error = 1; printf("[ERROR] missing band !\n");
  }
  if (opts->seg < 1) {
    error = 1; printf("[ERROR] missing segment !\n");
  }
  if (sett->nod < 1) {
    error = 1; printf("[ERROR] missing nod !\n");
  }
  if (sett->dt < 0.) {
    error = 1; printf("[ERROR] missing dt !\n");
  }

  if (error == 1) exit(EXIT_FAILURE);

  // test overlap and narrowdown
  // todo: NAV
  if (opts->narrowdown < 0.) {
    opts->narrowdown = (1. - opts->overlap)/2. ;
  }
  opts->narrowdown *= M_PI;
  
  sett->fpo = 10. + (1. - opts->overlap)*opts->band*(0.5/sett->dt);
  
}


/* Generate grid from the M matrix (grid.bin) */ 

void read_grid( Search_settings *sett, 
		Command_line_opts *opts ) {
     
     sett->M = (double *) calloc (16, sizeof (double));
     
     FILE *data;
     char filename[1024];
     int i, status;
     
     // In case when -usedet option is used for one detector
     // i.e. opts->usedet has a length of 2 (e.g. H1 or V1),
     // or data dir contains only one detector,
     // read grid.bin from this detector subdirectory 
     // (see detectors_settings() in settings.c for details)
     
     if(strlen(opts->usedet)==2 || sett->nifo==1 ) {
	  //sprintf (filename, "%s/%03d/%s/grid.bin", opts->dtaprefix, opts->ident, opts->usedet);
     	  sprintf (filename, "%s/%03d/%s/grid.bin", opts->indir, opts->seg, ifo[0].name);
     } else {
	  //sprintf (filename, "%s/%03d/grid.bin", opts->dtaprefix, opts->ident);
	  // det network
	  char dnet_str[MAX_DETECTORS*DETNAME_LENGTH]="";
          for(i=0; i<sett->nifo; i++)
	       strcat(dnet_str, ifo[i].name);

	  sprintf (filename, "%s/%03d/grids/grid_%03d_%04d_%sc.bin", opts->indir, opts->seg, opts->seg, opts->band, dnet_str);
     }
	  

     if ((data=fopen (filename, "r")) != NULL) {
	  printf("Using grid file %s\n", filename);
	  status = fread ((void *)&sett->fftpad, sizeof (int), 1, data);

	  printf("fftpad from the grid file: %d\n", sett->fftpad); 
	
	  // M: vector of 16 components consisting of 4 rows
	  // of 4x4 grid-generating matrix
	  status = fread ((void *)sett->M, sizeof (double), 16, data);
	  fclose (data);
     } else {
	  perror (filename);
	  exit(EXIT_FAILURE);
     }

} // end of read grid 


  /* Array initialization */ 

void init_arrays(
		 Search_settings *sett, 
		 Command_line_opts *opts,
		 Aux_arrays *aux_arr
		 ) {

  int i; 
  size_t status;

  // Allocates and initializes to zero the data, detector ephemeris
  // and the F-statistic arrays

  FILE *data;

  for(i=0; i<sett->nifo; i++) { 

    ifo[i].sig.xDat = (double *) calloc(sett->N, sizeof(double));
    
    // Extra array in case of software injections (to keep the original contents)
    //#mb  
    //ifo[i].sig.xDatorig = (double *) calloc(sett->N, sizeof(double));
    
    // Input time-domain data handling
    // 
    // The file name ifo[i].xdatname is constructed 
    // in settings.c, while looking for the detector 
    // subdirectories
    
    if((data = fopen(ifo[i].xdatname, "r")) != NULL) {
      status = fread((void *)(ifo[i].sig.xDat), 
		     sizeof(double), sett->N, data);
      fclose (data);
      
    } else {
      perror (ifo[i].xdatname);
      exit(EXIT_FAILURE); 
    }
    
    //#mb
    // Copy ifo[i].sig.xDat to ifo[i].sig.xDatorig for safekeeping 
    //memcpy(ifo[i].sig.xDatorig, ifo[i].sig.xDat, sett->N*sizeof(double));

    int j, Nzeros=0;
    // Checking for null values in the data
    for(j=0; j < sett->N; j++)
      if(!ifo[i].sig.xDat[j]) Nzeros++;
    
    ifo[i].sig.Nzeros = Nzeros; 
    
    // factor N/(N - Nzeros) to account for null values in the data
    ifo[i].sig.crf0 = (double)sett->N/(sett->N - ifo[i].sig.Nzeros);

    // Estimation of the variance for each detector 
    ifo[i].sig.sig2 = (ifo[i].sig.crf0)*var(ifo[i].sig.xDat, sett->N);

    ifo[i].sig.DetSSB = (double *) calloc(3*sett->N, sizeof(double));
    /* 
    const size_t array_bytes = 3*sett->N*sizeof(double);
    ifo[i].sig.DetSSB = NULL;
    if ( posix_memalign((void**)&ifo[i].sig.DetSSB, 32, array_bytes) ) exit (1);
    */

    // Ephemeris file handling
    char filename[562];
    sprintf (filename, "%s/%03d/%s/DetSSB.bin", 
        opts->indir, opts->seg, ifo[i].name);

    if((data = fopen(filename, "r")) != NULL) {
      // Detector position w.r.t Solar System Baricenter
      // for every datapoint
      status = fread((void *)(ifo[i].sig.DetSSB), 
               sizeof(double), 3*sett->N, data);

      // Deterministic phase defining the position of the Earth
      // in its diurnal motion at t=0 
      status = fread((void *)(&ifo[i].sig.phir), 
               sizeof(double), 1, data);

      // Earth's axis inclination to the ecliptic at t=0
      status = fread((void *)(&ifo[i].sig.epsm), 
               sizeof(double), 1, data);
      fclose (data);

      printf("Using %s as detector %s ephemerids...\n", filename, ifo[i].name);

    } else {
      perror (filename);
      return ;
    }

    // sincos 
    ifo[i].sig.sphir = sin(ifo[i].sig.phir);
    ifo[i].sig.cphir = cos(ifo[i].sig.phir);
    ifo[i].sig.sepsm = sin(ifo[i].sig.epsm);
    ifo[i].sig.cepsm = cos(ifo[i].sig.epsm);

    sett->sepsm = ifo[i].sig.sepsm; 
    sett->cepsm = ifo[i].sig.cepsm; 

    ifo[i].sig.xDatma = fftw_malloc(sett->N*sizeof(complex double));
    ifo[i].sig.xDatmb = fftw_malloc(sett->N*sizeof(complex double));

    ifo[i].sig.aa = (double *) calloc(sett->N, sizeof(double));
    ifo[i].sig.bb = (double *) calloc(sett->N, sizeof(double));

    ifo[i].sig.shft = (double *) calloc(sett->N, sizeof(double));
    ifo[i].sig.shftf = (double *) calloc(sett->N, sizeof(double));
 
  } // end loop for detectors 

  // Check if the ephemerids have the same epsm parameter
  for(i=1; i<sett->nifo; i++) {  
    if(!(ifo[i-1].sig.sepsm == ifo[i].sig.sepsm)) { 
      printf("The parameter epsm (DetSSB.bin) differs for detectors %s and %s. Aborting...\n", ifo[i-1].name, ifo[i].name); 
      exit(EXIT_FAILURE);
    } 

  } 

  // if all is well with epsm, take the first value 
  sett->sepsm = ifo[0].sig.sepsm;
  sett->cepsm = ifo[0].sig.cepsm;
  
  // Auxiliary arrays, Earth's rotation
  aux_arr->t2 = (double *) calloc(sett->N, sizeof (double));
  aux_arr->cosmodf = (double *) calloc(sett->N, sizeof (double));
  aux_arr->sinmodf = (double *) calloc(sett->N, sizeof (double));
  double omrt;

  for (i=0; i<sett->N; i++) {
    omrt = (sett->omr)*i;     // Earth angular velocity * dt * i
    aux_arr->t2[i] = sqr((double)i);
    aux_arr->cosmodf[i] = cos(omrt);
    aux_arr->sinmodf[i] = sin(omrt);
  }
  
} // end of init arrays 

void grid_sky_position(
		Search_settings *sett,
		Command_line_opts *opts) { 


} 
	

  /* Add signal to data   */ 

void add_signal(
		Search_settings *sett,
		Command_line_opts *opts,
		Aux_arrays *aux_arr,
		Search_range *s_range, 
    char *line) {

  int i, j, n, reffr; 
  double snr=0, sum = 0., h0=0, cof, d1; 
  double sigma_noise = 1.0;
  double be[2];
  double al1, al2, sinalt, cosalt, sindelt, cosdelt, phaseadd, shiftadd;
  double ph_o, psik, hoc, hop, iota, intrang[3]; 
  double nSource[3], sgnlo[7], sgnlol[4];
 
  char amporsnr[4];  

  do { 
    sscanf (line, "%s", amporsnr);
  } while ( strcmp(amporsnr, "amp")!=0 && strcmp(amporsnr, "snr")!=0 );


  /* check if signal file contains information on GW amplitude or SNR
    other parameters are grid ranges in f, s d and a, the reference segment number 
    frequency: sgnlo[0], spindown: sgnlo[1], declination: sgnlo[2], right ascension: sgnlo[3], 
    sgnlo[4]-sgnlo[6]: 3 variables related to orientation of the source (see sigen.c)
  */
  if(!strcmp(amporsnr, "amp")) { 
      sscanf (&line[3], "%le %d %d %d %d %d %le %le %le %le %le %le %le", 
        &h0, &s_range->gsize_f, &s_range->gsize_s, &s_range->gsize_m, &s_range->gsize_n, &reffr, 
        sgnlo, sgnlo+1, sgnlo+2, sgnlo+3, sgnlo+4, sgnlo+5, sgnlo+6); 
      printf("add_signal(): GW amplitude h0 is %le\n", h0); 
    } else if(!strcmp(amporsnr, "snr")) { 
      sscanf (&line[3], "%le %d %d %d %d %d %le %le %le %le %le %le %le", 
        &snr, &s_range->gsize_f, &s_range->gsize_s, &s_range->gsize_m, &s_range->gsize_n, &reffr, 
        sgnlo, sgnlo+1, sgnlo+2, sgnlo+3, sgnlo+4, sgnlo+5, sgnlo+6);  
      printf("add_signal(): GW (network) signal-to-noise ratio is %le\n", snr); 
    } else { 
      printf("Problem with the signal parameters (neither amp nor snr). Exiting...\n"); 
      exit(0); 
    } 

  printf("add_signal(): signal parameters: f=%le, s=%le, d=%le, a=%le, ph_o=%le, psik=%le, iota=%le\n", 
    sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3], sgnlo[4], sgnlo[5], sgnlo[6]);

  // Search-specific parametrization of freq. for software injections
  // i.e. shifting the freq. to current time segment 
  // sgnlo[0]: frequency, sgnlo[1]: frequency. derivative  
  sgnlo[0] += -2.*sgnlo[1]*(sett->N)*(reffr - opts->seg); 

  // Check if the signal is in band 
  if( sgnlo[0]<0 || sgnlo[0]>M_PI ) {
      printf("add_signal(): signal out of band f=%le s=%le\n", sgnlo[0], sgnlo[1]);
      return;
  }

  // Saving the injection frequency in s_range struct 
  s_range->freq_inj = sgnlo[0];

  cof = sett->oms + sgnlo[0]; 
  
  for(i=0; i<2; i++) sgnlol[i] = sgnlo[i]; 
  
  // Calculate the hemisphere and be vector 
  s_range->pmr[0] = ast2lin(sgnlo[3], sgnlo[2], C_EPSMA, be);

  sgnlol[2] = be[0]*cof;  
  sgnlol[3] = be[1]*cof; 
	 	
  // solving a linear system in order to translate 
  // spindown and sky positions (sgnlo parameters) 
  // into the grid coordinates 
  // Writes to: s_range->spndr[0], s_range->nr[0], s_range->mr[0]
  sda_to_grid(sett, s_range, sgnlol); 
 
  //#mb start 
  //
  // Grid positions
  al1 = s_range->nr[0]*sett->M[10] + s_range->mr[0]*sett->M[14];
  al2 = s_range->nr[0]*sett->M[11] + s_range->mr[0]*sett->M[15];
 
  // Change linear (grid) coordinates to real coordinates
  lin2ast(al1/sett->oms, al2/sett->oms, 
	  s_range->pmr[0], sett->sepsm, sett->cepsm,
	  &sinalt, &cosalt, &sindelt, &cosdelt);

  // calculate declination and right ascention
  // written in file as candidate signal sky positions
  sgnlo[2] = asin(sindelt);
  sgnlo[3] = fmod(atan2(sinalt, cosalt) + 2.*M_PI, 2.*M_PI);

  printf("Nearest grid point sky positions: DEC %lf, RA %lf\n", sgnlo[2], sgnlo[3]); 

  /*  
  // Calculate the hemisphere and be vector 
  s_range->pmr[0] = ast2lin(sgnlo[3], sgnlo[2], C_EPSMA, be);

  sgnlol[2] = be[0]*cof;  
  sgnlol[3] = be[1]*cof; 
	 	
  // solving a linear system for the second time in order to 
  // re-calculate the spindown position for the new sky positions 
  // which now match the grid point exactly 
  sda_to_grid(sett, s_range, sgnlol); 

  //#mb end 
  */ 

  // Define the grid ranges in which the signal will be looked for
  // +- grid ranges in each direction are defined by gsize_X
  s_range->spndr[1] = s_range->spndr[0] + s_range->gsize_s;
  s_range->spndr[0] -= s_range->gsize_s;
  s_range->nr[1] = s_range->nr[0] + s_range->gsize_n; 
  s_range->nr[0] -= s_range->gsize_n;
  s_range->mr[1] = s_range->mr[0] + s_range->gsize_m; 
  s_range->mr[0] -= s_range->gsize_m;
  // This variable sets the hemisphere number 
  s_range->pmr[1] = s_range->pmr[0];  

  // starting values for the loops for spindown, sky and hemispheres 
  s_range->sst = s_range->spndr[0];
  s_range->nst = s_range->nr[0];
  s_range->mst = s_range->mr[0];
  s_range->pst = s_range->pmr[0];

  printf("add_signal(): following grid range is used (spndr, nr, mr, pmr pairs)\n");
  printf("%f %f %f %f %f %f %f %f\n", \
   s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
   s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);

  // sgnlo[2]: declination, sgnlo[3]: right ascension 
  sindelt = sin(sgnlo[2]); 
  cosdelt = cos(sgnlo[2]); 
  sinalt = sin(sgnlo[3]);  
  cosalt = cos(sgnlo[3]); 

  // Intrinsic parameters of the signal: phase, polarization, inclination wobble angle iota 
  ph_o = sgnlo[4];    
  psik = sgnlo[5];    
  hoc =  cos(sgnlo[6]); // sgnlo[6] = iota = acos(hoc)   
  hop = (1. + hoc*hoc)/2.; 

  intrang[0] =  cos(2.*psik)*hop*cos(ph_o) - sin(2.*psik)*hoc*sin(ph_o) ;
  intrang[1] =  sin(2.*psik)*hop*cos(ph_o) + cos(2.*psik)*hoc*sin(ph_o) ;
  intrang[2] = -cos(2.*psik)*hop*sin(ph_o) - sin(2.*psik)*hoc*cos(ph_o) ;
  intrang[3] = -sin(2.*psik)*hop*sin(ph_o) + cos(2.*psik)*hoc*cos(ph_o) ;

  // To keep coherent phase between time segments  
  double phaseshift = sgnlo[0]*sett->N*(reffr - opts->seg)   
    + sgnlo[1]*pow(sett->N*(reffr - opts->seg), 2); 

  // Allocate arrays for added signal, for each detector 
  double **signadd = malloc((sett->nifo)*sizeof(double *));
  for(n=0; n<sett->nifo; n++)
    signadd[n] = malloc((sett->N)*sizeof(double));

  printf("%s %.8le %.8le %.8le %.8le %.8le %.8le %.8le\n", 
    opts->si_label, sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3], sgnlo[4], sgnlo[5], sgnlo[6]); 

  // Loop for each detector - sum calculations
  for(n=0; n<sett->nifo; n++) {
    
    modvir(sinalt, cosalt, sindelt, cosdelt,
	   sett->N, &ifo[n], aux_arr);

    nSource[0] = cosalt*cosdelt;
    nSource[1] = sinalt*cosdelt;
    nSource[2] = sindelt;
					
    for (i=0; i<sett->N; i++) {
      
      shiftadd = 0.; 					 
      for (j=0; j<3; j++)
      	shiftadd += nSource[j]*ifo[n].sig.DetSSB[i*3+j];		 
      
      // Phase 
      phaseadd = sgnlo[0]*i + sgnlo[1]*aux_arr->t2[i] 
        + (cof + 2.*sgnlo[1]*i)*shiftadd
        - phaseshift; 

      // The whole signal with 4 amplitudes and modulations 
      signadd[n][i] = intrang[0]*(ifo[n].sig.aa[i])*cos(phaseadd) 
                    + intrang[1]*(ifo[n].sig.aa[i])*sin(phaseadd) 
                    + intrang[2]*(ifo[n].sig.bb[i])*cos(phaseadd) 
                    + intrang[3]*(ifo[n].sig.bb[i])*sin(phaseadd);

      // Sum over signals
      sum += pow(signadd[n][i], 2.);
    
    } // data loop
   
  } // detector loop


  // Signal amplitude h0 from the snr 
  // (currently only makes sense for Gaussian noise with fixed sigma)
  if(snr)
    h0 = (snr*sigma_noise)/(sqrt(sum));

  // Loop for each detector - adding signal to data (point by point)  								
  for(n=0; n<sett->nifo; n++) {
    for (i=0; i<sett->N; i++) {

      // Adding the signal to the data vector 
      if(ifo[n].sig.xDat[i]) { 
//#mb 
//        ifo[n].sig.xDat[i] = ifo[n].sig.xDatorig[i] + h0*signadd[n][i];
        ifo[n].sig.xDat[i] += h0*signadd[n][i];

      } 

    } // data loop

  } // detector loop

  // printf("snr=%le h0=%le\n", snr, h0);

  // Free auxiliary 2d array 
  for(n=0; n<sett->nifo; n++) 
    free(signadd[n]);
  free(signadd);
 
} // add_signal()


void sda_to_grid(Search_settings *sett,
		Search_range *s_range, 
		double sgnlol[]) { 

  int i, s; 
  double *MM; 
  MM = (double *) calloc (16, sizeof (double));
  for(i=0; i<16; i++) MM[i] = sett->M[i] ;
  
  gsl_vector *x = gsl_vector_alloc (4);     
  
  gsl_matrix_view m = gsl_matrix_view_array (MM, 4, 4);
  gsl_matrix_transpose (&m.matrix) ; 
  gsl_vector_view b = gsl_vector_view_array (sgnlol, 4);
  gsl_permutation *p = gsl_permutation_alloc (4);
  
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
  
  s_range->spndr[0] = gsl_vector_get(x,1); 
  s_range->nr[0]    = gsl_vector_get(x,2);
  s_range->mr[0]    = gsl_vector_get(x,3);
  
  gsl_permutation_free (p);
  gsl_vector_free (x);
  free (MM);

} 




/* Search range */ 

void set_search_range(
		      Search_settings *sett, 
		      Command_line_opts *opts, 
		      Search_range *s_range) { 
  
  // Hemispheres (with respect to the ecliptic)
  if(opts->hemi) {
    s_range->pmr[0] = opts->hemi;
    s_range->pmr[1] = opts->hemi;

  } else {
    s_range->pmr[0] = 1;
    s_range->pmr[1] = 2;
  }

  // If the parameter range is invoked, the search is performed
  // within the range of grid parameters from an ascii file
  // ("-r range_file" from the command line)
  FILE *data;
  if (strlen (opts->range_file)) {

    if ((data = fopen(opts->range_file, "r")) != NULL) {

      int aqq = fscanf(data, "%f %f %f %f %f %f %f %f",
		       s_range->spndr, 1+s_range->spndr, 
		       s_range->nr, 1+s_range->nr, 
		       s_range->mr, 1+s_range->mr,
		       s_range->pmr, 1+s_range->pmr);
      
      /*
      //#mb commented-out for now - useful for tests 

      // the case when range file does not contain 8 integers
      // describing the grid ranges, but other values:
      // the pulsar position, frequency, and spindowns.
      if(range_status!=8) {

      rewind(data);
      range_status = fscanf (data, "%le %le %le %le %le %le %d",
      &pepoch, &alpha, &delta, &f0, &f1, &f2, &gsize);

      // GPS time of the first sample
      double gps1;
      sprintf (filename, "%s/%03d/starting_date", dtaprefix, ident);
      if ((data2 = fopen (filename, "r")) != NULL) {
      fscanf (data2, "%le", &gps1);
      fclose(data2);
      } else {
      perror (filename);
      return 1;
      }

      // Conversion of mjd to gps time
      double pepoch_gps = (pepoch - 44244)*86400 - 51.184;
      //			 gps1 = (gps1 - 44244)*86400 - 51.184;

      // Interpolation of ephemeris parameters to the starting time
      double *sgnlo;
      sgnlo = (double *) calloc (4, sizeof (double));

      double *be;
      be = (double *) calloc (2, sizeof (double));

      // ast2lin (auxi.c) returns the hemisphere number
      // and the vector be (used for sky position in linear coords.)
      pmr[0] = ast2lin(alpha, delta, epsm, be);

      sgnlo[0] = f0 + f1*(gps1 - pepoch_gps) + f2*pow(gps1 - pepoch_gps, 2)/2.;
      sgnlo[0] = 2*M_PI*2*sgnlo[0]*dt - oms;

      sgnlo[1] = f1 + f2*(gps1 - pepoch_gps);
      sgnlo[1] = M_PI*2*sgnlo[1]*dt*dt;

      sgnlo[2] = be[0]*(oms + sgnlo[0]);
      sgnlo[3] = be[1]*(oms + sgnlo[0]);

      // solving a linear system in order to translate
      // sky position, frequency and spindown (sgnlo parameters)
      // into the position in the grid

      gsl_vector *x = gsl_vector_alloc (4);
      int s;

      gsl_matrix_view m = gsl_matrix_view_array (M, 4, 4);
      gsl_matrix_transpose (&m.matrix) ;
      gsl_vector_view b = gsl_vector_view_array (sgnlo, 4);
      gsl_permutation *p = gsl_permutation_alloc (4);

      gsl_linalg_LU_decomp (&m.matrix, p, &s);
      gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

      spndr[0] = round(gsl_vector_get(x, 1));
      nr[0]		= round(gsl_vector_get(x, 2));
      mr[0]		= round(gsl_vector_get(x, 3));

      gsl_permutation_free (p);
      gsl_vector_free (x);
      free (be);
      free (sgnlo);


      // Warnings and infos
      if(hemi)

      printf("Warning: -h switch hemisphere choice (%d) may be altered\nby the choice of -r grid range...\n", hemi);

      // Define the grid range in which the signal will be looked for
      spndr[1] = spndr[0] + gsize ;
      spndr[0] -= gsize;
      nr[1] = nr[0] + gsize ;
      nr[0] -= gsize;
      mr[1] = mr[0] + gsize ;
      mr[0] -= gsize;
      pmr[1] = pmr[0];

      }
      */

      fclose (data);

    } else {
      perror (opts->range_file);
      exit(EXIT_FAILURE);
    }


  } else {

    // Establish the grid range in which the search will be performed
    // with the use of the M matrix from grid.bin
    gridr(
	  sett->M, 
	  s_range->spndr,
	  s_range->nr,
	  s_range->mr,
	  sett->oms,
	  sett->Smax);

    if (strlen(opts->dump_range_file)) {

      FILE *data;
      if ((data=fopen (opts->dump_range_file, "w")) != NULL) {
	fprintf(data, "%f %f\n%f %f\n%f %f\n%f %f\n",
		s_range->spndr[0], s_range->spndr[1],
		s_range->nr[0], s_range->nr[1],
		s_range->mr[0], s_range->mr[1],
		s_range->pmr[0], s_range->pmr[1] );
	
	printf("Wrote input data grid ranges to %s\n", opts->dump_range_file);
	fclose (data);
	exit(EXIT_SUCCESS);
	
      } else {
	
	printf("Can't open %s file for writing\n", opts->dump_range_file);
       	exit(EXIT_FAILURE);
	
      }
    }

  }

  printf("set_search_range() - the grid ranges are maximally this:\n"); 
  printf("(spndr, nr, mr, pmr pairs): %f %f %f %f %f %f %f %f\n",
	 s_range->spndr[0], s_range->spndr[1], s_range->nr[0], s_range->nr[1],
	 s_range->mr[0], s_range->mr[1], s_range->pmr[0], s_range->pmr[1]);

  printf("Smin: %le, -Smax: %le\n", sett->Smin, sett->Smax);

} // end of set search range 


  /* FFT Plans 
	 */

void plan_fftw(
	       Search_settings *sett,
	       Command_line_opts *opts,
	       FFTW_plans *plans,
	       FFTW_arrays *fftw_arr,
	       Aux_arrays *aux_arr) {

  char hostname[128], wfilename[512];
  FILE *wisdom;

  /* Imports a "wisdom file" containing information 
   * (previous tests) about how to optimally compute Fourier 
   * transforms on a given machine. If wisdom file is not present, 
   * it will be created after the test (measure) runs 
   * of the fft_plans are performed below 
   * (see http://www.fftw.org/fftw3_doc/Wisdom.html)
   */ 

#if defined(_OPENMP)
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#ifdef COMP_FLOAT
  fftwf_init_threads();
  FFTW_PRE(_plan_with_nthreads)(omp_get_max_threads());
#endif
#endif
  
  gethostname(hostname, 128);
  sprintf (wfilename, "wisdom-%s.dat", hostname);
  if((wisdom = fopen (wfilename, "r")) != NULL) {
    //fftw_import_wisdom_from_file(wisdom);
    //if (fftwf_import_wisdom_from_file(wisdom) == 0 ) exit(1);
    fclose (wisdom);
  }

  // moved to settings
  //sett->Ninterp = sett->interpftpad*sett->nfft; 
  //sett->nfftf = sett->fftpad*sett->nfft;

  // arrays xa,xb are used for in-place interpolation,
  // thus their length is max{fftpad*nfft, Ninterp}
  fftw_arr->arr_len = (sett->nfftf > sett->Ninterp
		       ? sett->nfftf : sett->Ninterp);

  fftw_arr->xa = (fftw_complex *)fftw_malloc(fftw_arr->arr_len*sizeof(fftw_complex));
  fftw_arr->xb = (fftw_complex *)fftw_malloc(fftw_arr->arr_len*sizeof(fftw_complex));

  printf("[fft plans] arr_len=%d   nfft=%d   fftpad=%d  Ninterp=%d\n",
	 fftw_arr->arr_len, sett->nfft, sett->fftpad, sett->Ninterp);

  printf("double plans...");
  // Change FFTW_MEASURE to FFTW_PATIENT for more optimized plan
  plans->pl_int = fftw_plan_dft_1d(sett->nfft, fftw_arr->xa, fftw_arr->xa,
				   FFTW_FORWARD, FFTW_ESTIMATE);
  plans->pl_inv = fftw_plan_dft_1d(sett->Ninterp, fftw_arr->xa, fftw_arr->xa,
				   FFTW_BACKWARD, FFTW_ESTIMATE);
  printf("done\n");

  printf("float plans...");
  fftw_arr->fxa = (FFTW_PRE(_complex) *)FFTW_PRE(_malloc)(sett->nfftf*sizeof(FFTW_PRE(_complex)));
  fftw_arr->fxb = (FFTW_PRE(_complex) *)FFTW_PRE(_malloc)(sett->nfftf*sizeof(FFTW_PRE(_complex)));

  plans->plan = FFTW_PRE(_plan_dft_1d)(sett->nfftf, fftw_arr->fxa, fftw_arr->fxa,
				       FFTW_FORWARD, FFTW_ESTIMATE);
  printf("done\n");

  // Generates a wisdom FFT file if there is none
  if((wisdom = fopen(wfilename, "r")) == NULL) {
    wisdom = fopen(wfilename, "w");
    fftw_export_wisdom_to_file(wisdom);
#ifdef COMP_FLOAT    
    fftwf_export_wisdom_to_file(wisdom);
#endif    
  }

  fclose (wisdom);

} // end of FFT plans 


  /* Checkpointing */

void read_checkpoints( Command_line_opts *opts, 
		       Search_range *s_range, 
		       int *FNum) {

  if(opts->checkp_flag) {

    // filename of checkpoint state file, depending on the hemisphere
    if(opts->hemi)
      sprintf(opts->state_file, "state_%03d_%04d%s_%d.dat",  
	      opts->seg, opts->band, opts->label, opts->hemi);
    else
      sprintf(opts->state_file, "state_%03d_%04d%s.dat", 
	      opts->seg, opts->band, opts->label);
    
    FILE *state;
    if((state = fopen(opts->state_file, "r")) != NULL) {

      // Scan the state file to get last recorded parameters
      if((fscanf(state, "%d %d %d %d %d", &s_range->pst, &s_range->mst,
		 &s_range->nst, &s_range->sst, FNum)) == EOF) {

	// This means that state file is empty (=end of the calculations)
	fprintf (stderr, "State file empty: nothing to do...\n");
	fclose (state);
	exit(EXIT_FAILURE);
	
      }

      fclose (state);

      // No state file - start from the beginning
    } else {
      s_range->pst = s_range->pmr[0];
      s_range->mst = s_range->mr[0];
      s_range->nst = s_range->nr[0];
      s_range->sst = s_range->spndr[0];
      *FNum = 0;
    } // if state

  } else {
    s_range->pst = s_range->pmr[0];
    s_range->mst = s_range->mr[0];
    s_range->nst = s_range->nr[0];
    s_range->sst = s_range->spndr[0];
    *FNum = 0;
  } // if checkp_flag
  
} // end reading checkpoints


  /* Cleanup & memory free  */

void cleanup(
	Search_settings *sett,
	Command_line_opts *opts,
	Search_range *s_range,
	FFTW_plans *plans,
	FFTW_arrays *fftw_arr,
	Aux_arrays *aux
	) {

  int i; 

  for(i=0; i<sett->nifo; i++) {
    free(ifo[i].sig.xDat);
//#mb    free(ifo[i].sig.xDatorig); 
    free(ifo[i].sig.xDatma);
    free(ifo[i].sig.xDatmb);
    free(ifo[i].sig.DetSSB);
    free(ifo[i].sig.aa);
    free(ifo[i].sig.bb);
    free(ifo[i].sig.shftf);
    free(ifo[i].sig.shft);
  } 
	
  free(aux->sinmodf);
  free(aux->cosmodf);
  free(aux->t2);
	
  fftw_free(fftw_arr->xa);
  fftw_free(fftw_arr->xb);  

  free(sett->M);
	
  FFTW_PRE(_destroy_plan)(plans->plan);
  fftw_destroy_plan(plans->pl_int);
  fftw_destroy_plan(plans->pl_inv);
  /*
  fftw_destroy_plan(plans->pl_int2);
  fftw_destroy_plan(plans->pl_inv2);
  fftw_destroy_plan(plans->plan2);
  */
  fftw_forget_wisdom();
  fftw_cleanup();


} // end of cleanup & memory free 



	/*	Command line options handling: coincidences  
	 */ 
	
void handle_opts_coinc( Search_settings *sett,
			Command_line_opts_coinc *opts,
			int argc,
			char* argv[]) {
	
  opts->wd=NULL;
	
  strcpy (opts->prefix, TOSTR(PREFIX));
  //strcpy (opts->dtaprefix, TOSTR(DTAPREFIX));

  // Initial value of the number of days is set to 0
  sett->nod = 0;

  // Default initial value of the data sampling time 
  sett->dt = 0.5;

  opts->help_flag=0;
  static int help_flag=0;

  // Default value of the minimal number of coincidences 
  opts->mincoin=3; 

  // Default value of the narrow-down parameter 
  opts->narrowdown=0.5; 

  // Default value of the cell shift: 0000 (no shifts)
  opts->shift=0;

  // Default value of the cell scaling: 1111 (no scaling)
  opts->scalef=4;
  opts->scales=4;
  opts->scaled=4;
  opts->scalea=4;

  sett->fpo = -1;
  
  // Default signal-to-noise threshold cutoff
  opts->snrcutoff=6;

  // Reading arguments 

  while (1) {
    static struct option long_options[] = {
      {"help", no_argument, &help_flag, 1},
      // Cell shifts  
      {"shift", required_argument, 0, 's'},
      // Cell scaling 
      {"scale", required_argument, 0, 'z'},
      // Reference frame number 
      {"refr", required_argument, 0, 'r'},
      // output directory
      {"output", required_argument, 0, 'o'},
      // input data directory
      {"data", required_argument, 0, 'd'},
      // frequency band number
      {"band", required_argument, 0, 'b'},
      // fpo value
      {"fpo", required_argument, 0, 'p'},
      // data sampling time 
      {"dt", required_argument, 0, 't'},
      // a file with input files (triggers + grids)
      {"infiles", required_argument, 0, 'i'},
      // Location of the reference frame grid
      {"refgrid", required_argument, 0, 'g'},
      // Minimal number of coincidences recorded in the output  
      {"mincoin", required_argument, 0, 'm'},
      // Narrow down the frequency band (+- the center of band) 
      {"narrowdown", required_argument, 0, 'n'},
      // Signal-to-noise threshold cutoff  
      {"snrcutoff", required_argument, 0, 'c'},
      // number of days in the time-domain segment 
      {"nod", required_argument, 0, 'y'},
      // band overlap
      {"overlap", required_argument, 0, 'v'},
      {0, 0, 0, 0}
    };

    if (help_flag) {

      printf("polgraw-allsky periodic GWs: search for concidences among candidates\n");
      printf("Usage: ./coincidences -[switch1] <value1> -[switch2] <value2> ...\n") ;
      printf("Switches are:\n\n");
      printf("-output       Output directory (default is ./coinc-results)\n");
      printf("-shift        Cell shifts in fsda directions (4 digit number, e.g. 0101, default 0000)\n");
      printf("-scale        Cell scaling in fsda directions (coma separated, e.g. 32,8,4,4, default 4,4,4,4)\n");
      printf("-refr         Reference frame number\n");
      printf("-band         Band number\n");
      printf("-fpo          Reference band frequency fpo value\n");
      printf("-dt           Data sampling time dt (default value: 0.5)\n");
      printf("-infile       File containing the list of trigger and grid files\n");      
      printf("-refgrid      Location of the reference frame grid\n");
      printf("-mincoin      Minimal number of coincidences recorded\n");
      printf("-narrowdown   Narrow-down the frequency band (range [0, 0.5] +- around center)\n");
      printf("-nod          Number of days\n");
      printf("-snrcutoff    Signal-to-noise threshold cutoff (default value: 6)\n");
      printf("-overlap      Band overlap, fpo=10+(1-overlap)*band/(dt*2) ; obligatory if band is used\n\n");
      
      printf("Also:\n\n");
      printf("--help		This help\n");

      exit (0);
    }

    int option_index = 0;
    int c = getopt_long_only (argc, argv, "p:o:s:z:r:t:g:m:n:c:y:b:v:i:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 'p':
	 sett->fpo = atof(optarg);
	 break;
    case 's': // Cell shifts 
	 opts->shift = atof(optarg);
	 break;
    case 'z': // Cell scaling
	 opts->scalef = atoi(strtok(optarg,","));
	 opts->scales = atoi(strtok(NULL,","));
	 opts->scaled = atoi(strtok(NULL,","));
	 opts->scalea = atoi(strtok(NULL,","));      
	 break;
    case 'r':
	 opts->refr = atoi(optarg);
	 break;
    case 'o':
	 strcpy(opts->prefix, optarg);
	 break;
    case 't':
	 sett->dt = atof(optarg);
	 break;
    case 'i':
	 strcpy(opts->infile, optarg);
	 break;
    case 'g':
	 strcpy(opts->refgrid, optarg);
	 break;
    case 'm':
	 opts->mincoin = atoi(optarg);
	 break;
    case 'n':
	 opts->narrowdown = atof(optarg);
	 break;
    case 'c':
	 opts->snrcutoff = atof(optarg);
	 break;
    case 'y':
	 sett->nod = atoi(optarg);
	 break;
    case 'b':
	 opts->band = atoi(optarg);
	 break;
    case 'v':
	 opts->overlap = atof(optarg);
	 break;
    case '?':
	 break;
    default:
	 break ;
    } /* switch c */
  } /* while 1 */

  // Putting the parameter in triggers' frequency range [0, pi] 
  opts->narrowdown *= M_PI; 

  // Check if sett->nod was set up, if not, exit
  if(!(sett->nod)) { 
    printf("Number of days not set... Exiting\n"); 
    exit(EXIT_FAILURE); 
  } 

  printf("Number of days is %d\n", sett->nod); 

  printf("The SNR threshold cutoff is %.12f, ", opts->snrcutoff); 
  printf("corresponding to F-statistic value of %.12f\n", 
    pow(opts->snrcutoff, 2)/2. + 2); 

  if(!(opts->band)) { 
       printf("Band is not set... Exiting\n"); 
       exit(EXIT_FAILURE);
  }
  if(!(opts->overlap)) { 
       printf("Band overlap is not set... Exiting\n"); 
       exit(EXIT_FAILURE);
  }
  
  printf("Band=%04d  Overlap=%f\n", opts->band, opts->overlap);

  // hemi must be decoded from filename
  opts->hemi = -1;
  
  // Starting band frequency:
  // fpo_val is optionally read from the command line
  // Its initial value is set to -1
  if (!(sett->fpo >= 0)) {
       if (opts->band > 0 && opts->overlap >=0.) {
	    sett->fpo = 10. + (1. - opts->overlap)*opts->band*(0.5/sett->dt);
       } else {
	    printf("Band AND overlap or fpo must be specified!\n");
	    exit(EXIT_FAILURE);
       }
  }
  printf("The reference frequency fpo is %f\n", sett->fpo);

  printf("Cell scaling factors are: %d %d %d %d\n", opts->scalef, opts->scales,
	 opts->scaled, opts->scalea);
  
} // end of command line options handling: coincidences  


void manage_grid_matrix( Search_settings *sett, char *gridfile ) {

  FILE *data;
  int status; 

  sett->M = (double *)calloc(16, sizeof (double));

  if ((data=fopen(gridfile, "r")) != NULL) {
       printf("Reading grid file: %s\n", gridfile);
       status = fread ((void *)&sett->fftpad, sizeof (int), 1, data);
       //printf("fftpad from the grid file: %d\n", sett->fftpad); 
       status = fread ((void *)sett->M, sizeof(double), 16, data);
       // We actually need the second (Fisher) matrix from grid.bin, 
       // hence the second fread: 
       status = fread ((void *)sett->M, sizeof(double), 16, data);
       fclose (data);
  } else {
       perror(gridfile);
       exit(EXIT_FAILURE);
  }


  // Calculating the eigenvectors and eigenvalues 
  gsl_matrix_view m = gsl_matrix_view_array(sett->M, 4, 4);

  gsl_vector *eval = gsl_vector_alloc(4);
  gsl_matrix *evec = gsl_matrix_alloc(4, 4);

  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(4); 
  gsl_eigen_symmv(&m.matrix, eval, evec, w);
  gsl_eigen_symmv_free(w);

  double eigval[4], eigvec[4][4]; 
  // Saving the results to the settings struct sett->vedva[][]
  int i, j;
  for(i=0; i<4; i++) { 
       eigval[i] = gsl_vector_get(eval, i); 
       gsl_vector_view evec_i = gsl_matrix_column(evec, i);
       
       for(j=0; j<4; j++)   
	    eigvec[j][i] = gsl_vector_get(&evec_i.vector, j);               
  }

  // This is an auxiliary matrix composed of the eigenvector 
  // columns multiplied by a matrix with sqrt(eigenvalues) on diagonal  
  for(i=0; i<4; i++) { 
       for(j=0; j<4; j++) { 
	    sett->vedva[i][j]  = eigvec[i][j]*sqrt(eigval[j]); 
	    //printf("%.12le ", sett->vedva[i][j]); 
       } 
//      printf("\n"); 
  }


/* 
  //#mb matrix generated in matlab, for tests 
  double _tmp[4][4] = { 
    {-2.8622034614137332e-001, -3.7566564762376159e-002, -4.4001551065376701e-012, -3.4516253934827171e-012}, 
    {-2.9591999145463371e-001, 3.6335210834374479e-002, 8.1252443441098394e-014, -6.8170555119669981e-014}, 
    {1.5497867603229576e-005, 1.9167007413107127e-006, 1.0599051611325639e-008, -5.0379548388381567e-008}, 
    {2.4410008440913992e-005, 3.2886518554938671e-006, -5.7338464150027107e-008, -9.3126913365595100e-009},
  };

  { int i,j; 
  for(i=0; i<4; i++) 
    for(j=0; j<4; j++) 
      sett->vedva[i][j]  = _tmp[i][j]; 
  }

  printf("\n"); 

  { int i, j; 
  for(i=0; i<4; i++) { 
    for(j=0; j<4; j++) {
        printf("%.12le ", sett->vedva[i][j]);
      }
      printf("\n"); 
  } 

 } 
*/ 

  gsl_vector_free (eval);
  gsl_matrix_free (evec);

} // end of manage grid matrix  
