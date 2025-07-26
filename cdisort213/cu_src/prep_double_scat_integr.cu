// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= prep_double_scat_integr () ==============*/

/*
       Prepares double scattering integration according to alternative
       Buras-Emde algorithm(201X).

                I N P U T   V A R I A B L E S

       nphase    number of angles for which original phase function
                     (ds->phase) is defined
       ntau
       nf        number of angular phase integration grid point
                     (zenith angle, theta)
       mu_phase  cos(theta) grid of phase function
       phas2     residual phase function

                O U T P U T   V A R I A B L E S

       mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
       neg_phas  index whether phas2 is negative
       norm_phas normalization factor for phase integration

                I N T E R N A L   V A R I A B L E S

       f_phas2_abs absolute value of integrated phase function
                      phas2
       f_phas2     cumulative integrated phase function phas2
       df          step length for calculating mu_eq

   Called by- c_new_intensity_correction
   Calls- c_dbl_vector, locate
 -------------------------------------------------------------------*/

void prep_double_scat_integr (int nphase, int ntau,
			      int           nf,
			      double       *mu_phase,
			      double       *phas2,
			      double       *mu_eq,
			      int          *neg_phas,
			      double       *norm_phas)
{
  int it=0, i=0, lu=0;
  double *f_phas2_abs=NULL;
  double f_phas2=0.0, df=0.0;

  f_phas2_abs = c_dbl_vector(0,nphase,"f_phas2_abs");

  for (lu=1; lu<=ntau; lu++) {

    /* calculate integral of |phas2| (f_phas2_abs) */

    F_PHAS2_ABS(1) = 0.0;
    for (it=2; it<=nphase; it++)
      F_PHAS2_ABS(it) = F_PHAS2_ABS(it-1) +
	( MUP(it) - MUP(it-1) ) * 0.5 *
	( fabs( PHAS2(it,lu) ) + fabs ( PHAS2(it-1,lu) ) );

    /* define mu grid which is equidistant in f_phas2_abs (mu_eq);
       find areas of negative phas2 (neg_phas);
       define normalization (norm_phas) */

    f_phas2 = 0.0;
    df = F_PHAS2_ABS(nphase) / (nf-1);
    MU_EQ(1,lu) = -1.0;

    if ( PHAS2(1,lu) > 0.0 )
      NEG_PHAS(1,lu) = FALSE;
    else
      NEG_PHAS(1,lu) = TRUE;

    it = 1;
    for (i=2; i<=nf-1; i++) {
      f_phas2 += df;

      while ( F_PHAS2_ABS(it+1) < f_phas2 )
	it++;

      MU_EQ(i,lu) = MUP(it)
	+ ( f_phas2 - F_PHAS2_ABS(it) ) /
	( F_PHAS2_ABS(it+1) - F_PHAS2_ABS(it) ) *
	( MUP(it+1) - MUP(it) );

      if ( PHAS2(it,lu) > 0.0 && PHAS2(it+1,lu) > 0.0 )
	NEG_PHAS(i,lu) = FALSE;
      else {
	if ( PHAS2(it,lu) < 0.0 && PHAS2(it+1,lu) < 0.0 )
	  NEG_PHAS(i,lu) = TRUE;
	else {
	  if ( PHAS2(it,lu) + ( f_phas2 - F_PHAS2_ABS(it) ) /
	       ( F_PHAS2_ABS(it+1) - F_PHAS2_ABS(it) ) *
	       ( PHAS2(it+1,lu) - PHAS2(it,lu) ) > 0.0 )
	    NEG_PHAS(i,lu) = FALSE;
	  else
	    NEG_PHAS(i,lu) = TRUE;
	}
      }

    } /* end for i<nf */

    MU_EQ(nf,lu) = 1.0;
    if ( PHAS2(nphase,lu) > 0.0 )
      NEG_PHAS(nf,lu) = FALSE;
    else
      NEG_PHAS(nf,lu) = TRUE;

    NORM_PHAS(lu) = F_PHAS2_ABS(nphase) / ( (nf-1) * M_PI );

  } /* end for lu<ntau */

  free(f_phas2_abs);
}

/*============================= end of prep_double_scat_integr() ========*/
