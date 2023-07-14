#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_legendre.h>
#include "util.h"
#include "input.h"

int main( int argc, char **argv )
{

  parameterBlock block;

  getInput( (const char **)argv, argc, &block );

  /* copy parameters */

  int do_planar = block.do_planar;

  double T = block.T;
  double r = block.r;
  double kc = block.kc;
  double chi = block.chi;
  double c0;
  if( !chi ) c0 = 0.0; else c0 = block.c0;

  double eta = block.eta;
  double Ap = block.Ap;
  double D0 = block.D0;

  double mode_max = block.mode_max;
  double rate = block.rate;
  int Nsteps = block.Nsteps;

  /* calculated parameters */
  
  double A = 4*M_PI*r*r;
  double rho0;
  if( !chi ) rho0 = 0.; else rho0 = 2.*chi/Ap;

  struct timeval tp;
  gettimeofday( &tp, NULL );
  int random_seed = tp.tv_usec;
  static const gsl_rng_type * rng_T;
  static gsl_rng * rng_r;
  gsl_rng_env_setup();
  rng_T = gsl_rng_default;
  rng_r = gsl_rng_alloc (rng_T);
  gsl_rng_set( rng_r, random_seed );

  /* Planar simulation */
  if( do_planar ){
    double q_max = mode_max;
    int P_max = ((q_max*r)-2);
    double *q = (double *)malloc( sizeof(double) * (P_max + 1));
    for( int i=0; i <= P_max; i++){
      q[i] = (i+2)/r;
    }

    double *Lambda = (double *)malloc( sizeof(double) * (P_max + 1));
    double *D = (double *)malloc( sizeof(double) * (P_max + 1));
    double *taum = (double *)malloc( sizeof(double) * (P_max + 1));
    double *taup = (double *)malloc( sizeof(double) * (P_max + 1));
    
    for( int i=0; i <= P_max; i++){
      Lambda[i] = 1/(4 * eta *q[i]);
    }
    for( int i=0; i <= P_max; i++){
      if( chi ) D[i]= rho0*D0*q[i]*q[i]; else D[i] = 0.0;
    }
    for( int i=0; i <= P_max; i++){
      taum[i] = 4. * eta / (kc * q[i]*q[i]*q[i]);
    }
    for( int i=0; i <= P_max; i++){
      taup[i] = 1. / (q[i]*q[i] * D0);
    }
    
    double *uq = (double *)calloc( (P_max+1)*2, sizeof(double));
    double *pq = (double *)calloc( (P_max+1)*2, sizeof(double));

    double *fuq = (double *)calloc( (P_max+1)*2, sizeof(double));
    double *fpq = (double *)calloc( (P_max+1)*2, sizeof(double));

    double *randm = (double *)calloc( (P_max+1)*2, sizeof(double));
    double *randp = (double *)calloc( (P_max+1)*2, sizeof(double));

    double dt = taum[P_max]/100;
    double tmax = Nsteps*rate / dt;
    printf("Time Parameters:\n");
    printf("/****************************/\n");
    printf("dt:      %e\n", dt);
    printf("Max time:      %f\n", tmax);
    printf("/****************************/\n");
    
    for( int t = 0; t < tmax; t++ )
      {
	for( int i = 0; i<=P_max; i++ )
	  {
	    randm[i * 2] = sqrt(T*dt*Lambda[i]*A)*gsl_ran_gaussian(rng_r, 1.0);
            randm[i * 2 + 1] = sqrt(T*dt*Lambda[i]*A)*gsl_ran_gaussian(rng_r, 1.0);
            randp[i * 2] = sqrt(q[i]*q[i]*D[i]*dt*A)*gsl_ran_gaussian(rng_r, 1.0);
	    randp[i * 2 + 1] = sqrt(q[i]*q[i]*D[i]*dt*A)*gsl_ran_gaussian(rng_r, 1.0);

	    fuq[i * 2] = - kc * q[i]*q[i]*q[i]*q[i] * uq[i * 2] + Ap * c0 * (kc/2.) * q[i]*q[i] * pq[i * 2];
            fuq[i * 2 + 1] = - kc * q[i]*q[i]*q[i]*q[i] * uq[i * 2 + 1] + Ap * c0 * (kc/2.) * q[i]*q[i] * pq[i * 2 + 1];

	    if( !chi )
	      {
	        fpq[i * 2] = 0;
                fpq[i * 2 + 1] = 0;
	      }
	    else if( chi )
	      { 
	        fpq[i * 2] = - pq[i * 2]*T*Ap/(2.*chi*(1-chi)) + Ap * c0 * (kc/2.) * q[i]*q[i] * uq[i * 2];
                fpq[i * 2 + 1] = - pq[i * 2 + 1]*T*Ap/(2.*chi*(1-chi)) + Ap * c0 * (kc/2.) * q[i]*q[i] * uq[i * 2 + 1];
	      }

	    uq[i * 2] += dt * Lambda[i] * fuq[i * 2] + randm[i * 2];
	    uq[i * 2 + 1] += dt * Lambda[i] * fuq[i * 2 + 1] + randm[i * 2 + 1];
	    pq[i * 2] += dt * D[i] * fpq[i * 2] / T + randp[i * 2];
	    pq[i * 2 + 1] += dt * D[i] * fpq[i * 2 + 1] / T + randp[i * 2 + 1];
	  }
	FILE *uqfile = fopen("uq.txt","w");
	if( t % int(rate/dt)  == 0)
	  {
	    printf("%d ", t);
	    for(int i = 0; i<=P_max; i++ )
	      {
		printf("%d %e %d %e ", 0, uq[i * 2], 0, uq[i * 2 + 1]);
	      }
            printf("\n");
	  }
	fclose(uqfile);
	FILE *pqfile = fopen("pq.txt","w");
        if( t % int(rate/dt)  == 0)
          {
            printf("%d ", t);
            for(int i = 0; i<=P_max; i++ )
              {
                printf("%d %e %d %e ", 0, pq[i * 2], 0, pq[i * 2 + 1]);
              }
            printf("\n");
          }
	fclose(pqfile);
      }
    free(q);
    free(Lambda);
    free(D);
    free(taum);
    free(taup);
    free(uq);
    free(pq);
    free(randm);
    free(randp);
    free(fuq);
    free(fpq);
  }

  /* Vesicle Simulation */
  else if( !do_planar ){
    int l_max = mode_max - 2;
    int m_max = 2 * mode_max + 1;
    int *l = (int*)malloc(sizeof(int)*(l_max + 1));
    for( int i=0; i <= l_max; i++ ){
      l[i] = i+2;
    }

    int *m = (int*)malloc(sizeof(int)*(m_max + 1));
    for( int i=0; i <= m_max; i++ ){
      m[i] = ceil(i/2);
    }

    double *Lambda = (double *)malloc( sizeof(double) * (l_max + 1));
    double *D = (double *)malloc( sizeof(double) * (l_max + 1));
    double *taum = (double *)malloc( sizeof(double) * (l_max + 1));
    double *taup = (double *)malloc( sizeof(double) * (l_max + 1));
  
    for( int i=0; i <= l_max; i++ ){
      Lambda[i] = (1/(eta*r*r*r)) * (l[i]*(l[i]+1.)/((2.*l[i]+1.)*(2.*l[i]*l[i]+2.*l[i]-1.)));
    }
    for( int i=0; i <= l_max; i++ ){
      if( chi ) D[i] = (2*l[i]*(l[i]+1)*chi*(1-chi)*D0)/(Ap*r*r*r*r); else D[i] = 0.0;
    }
    for( int i=0; i <= l_max; i++ ){
      taum[i] = (1/(kc*Lambda[i]*l[i]*(l[i]+1)*(l[i]+2)*(l[i]-1)));
    }
    for( int i=0; i <= l_max; i++ ){
      taup[i] = r*r / (l[i]*(l[i]+1) * D0);
    }

    double *ulm = (double *)calloc( (l_max + 1)*(m_max + 1), sizeof(double));
    double *rholm = (double *)calloc( (l_max + 1)*(m_max + 1), sizeof(double));
    
    double *plm = (double *)malloc( sizeof(double) *  (l_max + 1)*(m_max + 1));
    for( int i = 0; i <= l_max; i++ )
      {
	for(int j = 1; j <= m_max; j++ )
	  {
	    if( m[j] > l[i]){
	      plm[m_max*i + j] = 0.;
	    }
	    else{
	      plm[i * m_max + j] = gsl_sf_legendre_sphPlm(l[i], m[j], cos(0.5*M_PI));
	    }
	  }
      }

    double *randm = (double *)calloc( (l_max + 1)*(m_max + 1), sizeof(double));
    double *randp = (double *)calloc( (l_max + 1)*(m_max + 1), sizeof(double));

    double *fulm = (double *)calloc( (l_max + 1)*(m_max + 1), sizeof(double));
    double *frholm = (double *)calloc( (l_max + 1)*(m_max + 1), sizeof(double));

    double dt = taum[l_max]/100;
    int tmax = Nsteps*rate / dt;
    printf("Time Parameters:\n");
    printf("/****************************/\n");
    printf("dt:      %e\n", dt);
    printf("Max time:      %d\n", tmax);
    printf("/****************************/\n");

    for( int t = 0; t < tmax; t++ )
      {
	double *vq = (double *)calloc((m_max + 1), sizeof(double));
	double *rhoq = (double *)calloc((m_max + 1), sizeof(double));
	for( int i = 0; i <= l_max; i++ )
	  {
	    for(int j = 1; j <= m_max; j++ )
	      {
		randm[i * m_max + j] = sqrt(2*T*dt*Lambda[i])*gsl_ran_gaussian(rng_r, 1.0);
		randp[i * m_max + j] = sqrt(2*D[i]*dt)*gsl_ran_gaussian(rng_r, 1.0);

		fulm[i * m_max + j] = - Lambda[i] * (kc * (l[i]+2) * (l[i]-1) * (l[i] *(l[i]+1)) * ulm[i * m_max + j] + kc*c0*Ap*r*(-2+l[i]+l[i]*l[i])*rholm[i * m_max + j]/2.);
		if( !chi )
		  {
		    frholm[i * m_max + j] = 0;
		  }
		else if( chi )
		  {
		    frholm[i * m_max + j] = - (D[i]/T) * ((T*Ap*r*r)/(2*chi*(1 - chi)) * rholm[i * m_max + j] + kc*c0*Ap*r*(-2+l[i]+l[i]*l[i])*ulm[i * m_max + j]/2.);
		  }
		ulm[i * m_max + j] += dt * fulm[i * m_max + j] + randm[i * m_max + j];
		rholm[i * m_max + j] += dt * frholm[i * m_max + j] + randp[i * m_max + j];

		vq[j] += ulm[i * m_max + j]*plm[i * m_max + j];
		rhoq[j] += rholm[i * m_max + j]*plm[i * m_max + j];
	      }
	  }
	FILE *uqfile = fopen("uq.txt","w");
	if( t % int(rate/dt)  == 0)
          {
            printf("%d ", t);
	    for(int i = 4; i <= m_max; i=i+2)
	      {
		printf("%d %e %d %e ", 0, vq[i], 0, vq[i+1]);
	      }
	    printf("\n");
	  }
	fclose(uqfile);
	FILE *pqfile = fopen("pq.txt","w");
        if( t % int(rate/dt)  == 0)
          {
            printf("%d ", t);
            for(int i = 4; i <= m_max; i=i+2)
              {
                printf("%d %e %d %e ", 0, rhoq[i], 0, rhoq[i+1]);
              }
            printf("\n");
          }
	fclose(pqfile);
      }
    free(l);
    free(m);
    free(Lambda);
    free(D);
    free(taum);
    free(taup);
    free(ulm);
    free(rholm);
    free(randm);
    free(randp);
    free(fulm);
    free(frholm);
  }
  gsl_rng_free (rng_r);
}

