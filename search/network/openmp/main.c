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
#include <time.h>
#include <dirent.h>
#include <signal.h>

#include "auxi.h"
#include "struct.h"
#include "settings.h"
#include "jobcore.h"
#include "init.h"

#ifndef CODEVER
#define CODEVER unknown
#endif


Detector_settings ifo[MAX_DETECTORS];
volatile sig_atomic_t save_state = 0;


int main (int argc, char* argv[]) {

  Command_line_opts opts;
  Search_settings sett;
  Search_range s_range; 
  Aux_arrays aux_arr;
  int i; 

#define QUOTE(x) #x
#define STR(macro) QUOTE(macro)
#define CVSTR STR(CODEVER)

  printf("Code version : " CVSTR "\n");
  if (signal(SIGUSR1, sig_handler) != SIG_ERR &&
      signal(SIGTERM, sig_handler) != SIG_ERR )
    printf("State saved on SIGTERM or SIGUSR1\n");    

  // Command line options 
  read_ini_file(&sett, &opts, argc, argv);    

  // Directed search check 
  if (strlen(opts.ra) > 0 && strlen(opts.dec) > 0) {
    opts.is_directed = 1; 

    fprintf(stderr, "\nRA and DEC are both specified in the ini file - this is now a directed search.\n");
    opts.ra_val = atof(opts.ra);
    opts.dec_val = atof(opts.dec);

  } else { 
    opts.is_directed = 0; 
  } 

  // testing printout 
  if (opts.is_directed) {
    printf("%lf %lf\n", opts.ra_val, opts.dec_val);
  }  
  exit(0);

  // Output data handling
  struct stat buffer;

  if (stat(opts.outdir, &buffer) == -1) {
    if (errno == ENOENT) {
      // Output directory apparently does not exist, try to create one
      if(mkdir(opts.outdir, S_IRWXU | S_IRGRP | S_IXGRP 
          | S_IROTH	| S_IXOTH) == -1) {
	      perror (opts.outdir);
	      return 1;
      }
    } else { // can't access output directory
      perror (opts.outdir);
      return 1;
    }
  }
 
  // Detector network settings
  detectors_settings(&sett, &opts); 

  // Grid data 
  read_grid(&sett, &opts);	

  // Search settings
  search_settings(&sett); 

  // Array initialization, reading in the input data and the ephemerids 
  init_arrays(&sett, &opts, &aux_arr);

  // Narrowing-down the band (excluding the edges 
  // according to the opts.narrowdown parameter)
  // adds two lines
  if(opts.narrowdown < 0.5*M_PI) narrow_down_band(&sett, &opts);

/* #mb commented out for now   
  // Reading veto lines data from external files
  printf("Reading veto files...\n");
  read_lines(&sett, &opts);
  if (opts.gen_vlines_flag) exit(EXIT_SUCCESS);
*/ 

  // Amplitude modulation functions for each detector  
  for(i=0; i<sett.nifo; i++)   
    rogcvir(&ifo[i]); 

  // FFT plans 
  FFTW_plans fftw_plans;
  FFTW_arrays fftw_arr;
  plan_fftw(&sett, &opts, &fftw_plans, &fftw_arr, &aux_arr);

  // If addsig switch is used, signal(s) parameters are read from file, 
  // and search around signal's position (+- gsize, defined in file) is performed 
  if(strlen(opts.addsig)) { 

    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(opts.addsig, "r");
    if (fp == NULL) exit(EXIT_FAILURE);

    // signal counter (line number in file) 
    int sid = 1;     

    while ((read = getline(&line, &len, fp)) != -1) {

        sprintf(opts.si_label, "sig%03d", sid);

        // search settings, adding signal to data 
        // and corresponding grid range around the injection   
        add_signal(&sett, &opts, &aux_arr, &s_range, line);

        int Fnum=0;	// candidate signal number

        // main search job
        search(&sett, &opts, &s_range, 
            &fftw_plans, &fftw_arr, &aux_arr, &Fnum);

        sid++; 

    } 

  // normal run 
  } else {  

    // Set search range (either maximal range or from an optional range file)  
    set_search_range(&sett, &opts, &s_range);

    // Checkpointing
    int Fnum=0;	// candidate signal number
    read_checkpoints(&opts, &s_range, &Fnum);

    // main search job
    search(&sett, &opts, &s_range, 
         &fftw_plans, &fftw_arr, &aux_arr, &Fnum);
  } 

  // state file is emptied and closed in jobcore to mark successful end
  // do not remove it here
  /*
  FILE *state;
  if(opts.checkp_flag) {
    remove(opts.qname);
  }
  */
  
  // Cleanup & memory free 
  cleanup(&sett, &opts, &s_range, &fftw_plans, &fftw_arr, &aux_arr);

  return 0; 
	
}


// signal handler to save state and exit before end
static void sig_handler(int signo)
{
  if (signo == SIGTERM || signo == SIGUSR1) save_state = 1;
}
