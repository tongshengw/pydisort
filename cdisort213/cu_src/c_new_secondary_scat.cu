// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_new_secondary_scat() ==================*/

/*
   Calculates secondary scattered intensity, new method (see BDE)

                I N P U T   V A R I A B L E S

        ds        Disort state variables
        iu        index of user polar angle
        lu        index of user level
	it	  index where ctheta contained in mu grid of exact
	             phase function
        ctheta    cosine of scattering angle
        flyr      separated fraction f in Delta-M method
        layru     index of utau in multi-layered system
        tauc      cumulative optical depth at computational layers
        nf        number of angular phase integration grid point
                     (zenith angle, theta)
        phas2     residual phase function
        mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
        neg_phas  index whether phas2 is negative
        norm_phas normalization factor for phase integration

                I N T E R N A L   V A R I A B L E S

        pspike  2*p"-p"*p", where p" is the residual phase function
        pspike1 2*p", where p" is the residual phase function
        pspike2 p"*p", where p" is the residual phase function
        wbar    mean value of single scattering albedo
        fbar    mean value of separated fraction f
        dtau    layer optical depth
        stau    sum of layer optical depths between top of atmopshere and layer layru
	umu0p
        nphase  number of angles for which original phase function
                   (ds->phase) is defined

   Called by- c_new_intensity_correction
   Calls- calc_phase_squared, c_xi_func
 -------------------------------------------------------------------*/

double c_new_secondary_scat(disort_state *ds,
			    int           iu,
			    int           lu,
			    int           it,
			    double        ctheta,
			    double       *flyr,
			    int           layru,
			    double       *tauc,
			    int           nf,
			    double       *phas2,
			    double       *mu_eq,
			    int          *neg_phas,
			    double        norm_phas)
{
  register int
    lyr;
  const double
    tiny = 1.e-4;
  double
    dtau,fbar,pspike,
    stau,umu0p,wbar;
  int nphase=ds->nphase;

  double pspike1=0.0, pspike2=0.0;

  /*
   * Calculate vertically averaged value of single scattering albedo and separated
   * fraction f, eq. STWL (A.15)
   */
  dtau = UTAU(lu)-TAUC(layru-1);
  wbar = SSALB(layru)*dtau;
  fbar = FLYR(layru)*wbar;
  stau = dtau;
  for (lyr = 1; lyr <= layru-1; lyr++) {
    wbar += DTAUC(lyr)*SSALB(lyr);
    fbar += DTAUC(lyr)*SSALB(lyr)*FLYR(lyr);
    stau += DTAUC(lyr);
  }

  if (wbar <= tiny || fbar <= tiny || stau <= tiny || ds->bc.fbeam <= tiny) {
    return 0.;
  }

  fbar /= wbar;
  wbar /= stau;

  /* Calculate pspike1=P" */

  pspike1 = PHAS2(it,lu) + ( ctheta - ds->MUP(it) ) /
    ( ds->MUP(it+1) - ds->MUP(it) ) * ( PHAS2(it+1,lu) - PHAS2(it,lu) );

  pspike2 = calc_phase_squared (ds->nphase, lu, ctheta, nf,
				ds->mu_phase, phas2, mu_eq, neg_phas,
				norm_phas);

  pspike = 2.*pspike1 - pspike2;

  umu0p = ds->bc.umu0/(1.-fbar*wbar);

  /*
   * Calculate IMS correction term, eq. STWL (A.13)
   */
  return ds->bc.fbeam/(4.*M_PI)*SQR(fbar*wbar)/(1.-fbar*wbar)*pspike*c_xi_func(-UMU(iu),umu0p,UTAU(lu));
}

/*============================= end of c_new_secondary_scat() ===========*/
