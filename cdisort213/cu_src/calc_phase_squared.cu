// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= calc_phase_squared() ====================*/

/*
   Calculates squared phase function (see BDE)

                I N P U T   V A R I A B L E S

        nphase  number of angles for which original phase function
                   (ds->phase) is defined
        lu        index of user level
        ctheta  cosine of scattering angle
        nf        number of angular phase integration grid point
                     (zenith angle, theta)
        mu_phase  cos(theta) grid of phase function
        phas2     residual phase function
        mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
        neg_phas  index whether phas2 is negative
        norm_phas normalization factor for phase integration

                I N T E R N A L   V A R I A B L E S

        pspike2  p"*p", where p" is the residual phase function; return value
	mu1arr
	stheta   corresponding sin of ctheta
	smueq    corresponding sin of mu_eq
	phint    phase function integrated over phi
	scr

   Called by- c_new_secondary_scat
 -------------------------------------------------------------------*/

double calc_phase_squared (int           nphase,
			   int           lu,
			   double        ctheta,
			   int           nf,
			   double       *mu_phase,
			   double       *phas2,
			   double       *mu_eq,
			   int          *neg_phas,
			   double        norm_phas)
{
  int j=0, k=0, it=0;

  double pspike2=0.0, stheta=0.0;
  double smueq=0.0, phint=0.0;

  double mumin=0.0, mumax=0.0;
  int imin=0, imax=0;
  double D=0.0, C=0.0, Dp=0.0, Cp=0.0;
  int cutting=FALSE;

  stheta = sqrt( 1.0 - ctheta * ctheta );

  /* calculate pspike2 */

  /* Note: MU_EQ(j.lu) is mu_1; ctheta is mu; MUP(k) is mu_i in BDE(201X) */


  for (j=1;j<=nf;j++) {

    /* special case: second scattering angle does not depend on
       azimuth of first scattering angle */
    if (ctheta==1.0 || MU_EQ(j,lu)==1.0) {
      it = locate_disort ( mu_phase, nphase, MU_EQ(j,lu)*ctheta ) + 1;
      phint = M_PI * ( PHAS2(it,lu)
		       + ( MU_EQ(j,lu)*ctheta - MUP(it) )
		       / ( MUP (it+1) - MUP(it) )
		       * ( PHAS2(it+1,lu) - PHAS2(it,lu) ) );
      if (ctheta==1.0)
	phint /= 2.0;
    }
    else {
      phint = 0.0;

      smueq = sqrt ( 1. - MU_EQ(j,lu)*MU_EQ(j,lu) );

      /* locate integration borders */
      mumin = ctheta *  MU_EQ(j,lu) - stheta * smueq;
      mumax = ctheta *  MU_EQ(j,lu) + stheta * smueq;

      /* cut where mu_1 = mu_2 */
      if (MU_EQ(j,lu) < mumax) {
	mumax = MU_EQ(j,lu);
	cutting=TRUE;
      }
      else
	cutting=FALSE;

      if (mumin<mumax) {
	imin = locate_disort ( mu_phase, nphase, mumin)+1;
	imax = locate_disort ( mu_phase, nphase, mumax)+1;

	k=imin;
	/* assuming SPF is linear in mu */
	D = ( PHAS2(k+1,lu) - PHAS2(k,lu) ) / ( MUP(k+1) - MUP(k) );
	C = PHAS2(k,lu) - MUP(k) * D;

	phint +=  ( D * ctheta * MU_EQ(j,lu) + C ) * M_PI / 2.0;

	for (k=imin+1;k<=imax;k++) {

	  Dp = ( PHAS2(k+1,lu) - PHAS2(k,lu) ) / ( MUP(k+1) - MUP(k) );
	  Cp = PHAS2(k,lu) - MUP(k) * Dp;

	  phint +=
	    ( Dp - D ) * sqrt ( 1.0 - ctheta * ctheta
				- MU_EQ(j,lu) * MU_EQ(j,lu)
				+ 2.0 * ctheta * MU_EQ(j,lu) * MUP(k)
				- MUP(k) * MUP(k) )
	    + ( ( Dp - D )* ctheta * MU_EQ(j,lu) + Cp - C ) *
	    asin ( ( ctheta * MU_EQ(j,lu) - MUP(k) )
		   / ( smueq * stheta ) );

	  D=Dp;
	  C=Cp;
	}

	if (cutting==TRUE)
	  phint += - D * sqrt ( 1.0 - ctheta * ctheta
			      + 2.0 * MU_EQ(j,lu) * MU_EQ(j,lu) *
			      ( ctheta - 1.0 ) )
	    - ( D * ctheta * MU_EQ(j,lu) + C ) *
	    asin ( ( ctheta - 1.0 ) * MU_EQ(j,lu)
		   / ( smueq * stheta ) );
	else
	  phint += ( D * ctheta * MU_EQ(j,lu) + C ) * M_PI / 2.0;
      }
    }

    if (j==1 || j==nf) {
      if ( NEG_PHAS(j,lu) == TRUE )
	pspike2 = pspike2 - 0.5 * phint;
      else
	pspike2 = pspike2 + 0.5 * phint;
    }
    else {
      if ( NEG_PHAS(j,lu) == TRUE )
	pspike2 = pspike2 - phint;
      else
	pspike2 = pspike2 + phint;
    }

  }

  pspike2 *= norm_phas;

  return pspike2;
}

/*============================= end of calc_phase_squared() =============*/
