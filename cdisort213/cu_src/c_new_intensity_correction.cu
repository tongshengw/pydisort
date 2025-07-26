// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_new_intensity_correction() ============*/

/*
       Corrects intensity field by using alternative Buras-Emde algorithm
       (201X).

                I N P U T   V A R I A B L E S

       ds      Disort state variables
       dither  small multiple of machine precision
       flyr    separated fraction in delta-M method
       layru   index of UTAU in multi-layered system
       lyrcut  logical flag for truncation of computational layer
       ncut    total number of computational layers considered
       oprim   delta-M-scaled single-scatter albedo
       phirad  azimuthal angles in radians
       tauc    optical thickness at computational levels
       taucpr  delta-M-scaled optical thickness
       utaupr  delta-M-scaled version of UTAU

                O U T P U T   V A R I A B L E S

       out->UU  corrected intensity field; UU(IU,LU,J)
                 iu=1,ds->numu; lu=1,ds->ntau; j=1,ds->nphi

                I N T E R N A L   V A R I A B L E S

       ctheta    cosine of scattering angle
       dtheta    angle (degrees) to define aureole region as
                      direction of beam source +/- DTHETA
       phasa     actual (exact) phase function
       phasm     delta-M-scaled phase function
       phast     phase function used in TMS correction; actual phase
                      function divided by (1-FLYR*SSALB)
       pl        ordinary Legendre polynomial of degree l, P-sub-l
       plm1      ordinary Legendre polynomial of degree l-1, P-sub-(l-1)
       plm2      ordinary Legendre polynomial of degree l-2, P-sub-(l-2)
       theta0    incident zenith angle (degrees)
       thetap    emergent angle (degrees)
       ussndm    single-scattered intensity computed by using exact
                     phase function and scaled optical depth
                     (first term in STWL(68a))
       ussp      single-scattered intensity from delta-M method
                     (second term in STWL(68a))
       duims     intensity correction term from IMS method
                     (delta-I-sub-IMS in STWL(A.19))
       nf        number of angular phase integration grid point
                     (zenith angle, theta)
       np        number of angular phase integration grid point
                     (azimuth angle, phi)
       nphase    number of angles for which original phase function
                     (ds->phase) is defined
       mu_eq     cos(theta) phase integration grid points,
                     equidistant in abs(f_phas2)
       norm_phas normalization factor for phase integration
       norm      normalization factor for preparation of phas2
       neg_phas  index whether phas2 is negative
       phas2     residual phase function
       phasr     delta-M scaled phase function
       f_phas2   cumulative integrated phase function phas2
       fbar      mean value of separated fraction f

   Called by- c_disort
   Calls- c_single_scat, c__new_secondary_scat,
          prep_double_scat_integr, c_dbl_vector
 -------------------------------------------------------------------*/

void c_new_intensity_correction(disort_state  *ds,
				disort_output *out,
				double         dither,
				double        *flyr,
				int           *layru,
				int            lyrcut,
				int            ncut,
				double        *oprim,
				double        *phasa,
				double        *phast,
				double        *phasm,
				double        *phirad,
				double        *tauc,
				double        *taucpr,
				double        *utaupr)
{
  register int
    iu,jp,k,lc,ltau,lu;
  double
    ctheta,dtheta,duims,pl,plm1,plm2,
    theta0=0,thetap=0,ussndm,ussp;

  const int
    nf = 100;
  const double
    tiny = 1e-4;
  int it=0, lyr=0;
  int nphase=ds->nphase;

  double *mu_eq=NULL, *norm_phas=NULL, norm=0.0;
  int *neg_phas=NULL;

  double *phas2=NULL, *phasr=NULL;
  double f_phas2=0.0;
  double fbar=0.0;
  int need_secondary_scattering=0;

  dtheta = 10.;

  /* beginning of BDE stuff */

  /* check whether secondary scattering is performed at all */
  for (iu = 1; iu <= ds->numu; iu++) {
    if (UMU(iu) < 0.) {
      /*
       * Calculate zenith angles of incident and emerging directions
       */
      theta0 = acos(-ds->bc.umu0)/DEG;
      thetap = acos(UMU(iu))/DEG;
      if (fabs(theta0-thetap) <= dtheta) {
	need_secondary_scattering=TRUE;
	break;
      }
    }
  }

  if (need_secondary_scattering==TRUE) {
    /* Initialization of new PSPIKE.                                      */

    mu_eq  = c_dbl_vector(0,nf*ds->ntau-1,"mu_eq");
    norm_phas = c_dbl_vector(0,ds->ntau-1,"norm_phas");
    neg_phas  = c_int_vector(0,nf*ds->ntau-1,"neg_phas");
    phas2 = c_dbl_vector(0,ds->nphase*ds->ntau-1,"phas2");
    phasr = c_dbl_vector(0,ds->nlyr-1,"phasr");

    /* Calculate delta-scaled phase function (phasr) */

    for (it=1; it<=ds->nphase; it++) {

      ctheta = ds->MUP(it);

      for (lc=1; lc<=ds->nlyr; lc++)
	PHASR(lc) = 1.0 - FLYR(lc);

      plm1 = 1.0;
      plm2 = 0.0;

      for (k=1; k<=ds->nstr-1; k++) {

	/* ** Calculate Legendre polynomial of */
	/* ** P-sub-l by upward recurrence     */

	pl = ( (2*k-1) * ctheta * plm1 - (k-1) * plm2 ) / k;
	plm2 = plm1;
	plm1 = pl;

	for (lc=1; lc<=ds->nlyr; lc++)
	  PHASR(lc) += (2*k+1) * pl * ( PMOM(k,lc) - FLYR(lc) );

      }

      /* calculate difference between original and delta-scaled phase
	 functions (phas2) */

      for (lu=1; lu<=ds->ntau; lu++) {

	PHAS2(it,lu) = 0.0;

	/* this could be optimized */
	for (lyr=1; lyr<=LAYRU(lu)-1; lyr++)
	  PHAS2(it,lu) += ( DSPHASE(it,lyr) - PHASR(lyr) ) *
	    SSALB(lyr) * DTAUC(lyr);

	lyr = LAYRU(lu);
	PHAS2(it,lu) += ( DSPHASE(it,lyr) - PHASR(lyr) ) *
	  SSALB(lyr) * ( UTAU(lu) - TAUC(lyr-1) );

      }

    } /* end for it<nphas */

    /* normalize by 1/(ssa*beta*f) */

    for (lu=1; lu<=ds->ntau; lu++) {

      lyr = LAYRU(lu);
      fbar = FLYR(lyr) * SSALB(lyr) * ( UTAU(lu) - TAUC(lyr-1) );

      for (lyr=1; lyr<=LAYRU(lu)-1; lyr++)
	fbar += SSALB(lyr) * DTAUC(lyr) * FLYR(lyr);

      if ( fbar <= tiny || ds->bc.fbeam <= tiny )
	for (it=1; it<=ds->nphase; it++)
	  PHAS2(it,lu) = 0.0;
      else {
	fbar = 1. / fbar;
	for (it=1; it<=ds->nphase; it++)
	  PHAS2(it,lu) *= fbar;
      }

      /* normalize phas2 to 2.0 */

      f_phas2 = 0.0;
      for (it=2; it<=ds->nphase; it++)
	f_phas2 +=
	  ( ds->MUP(it) - ds->MUP(it-1) ) * 0.5 *
	  ( PHAS2(it,lu) + PHAS2(it-1,lu) );

      if (f_phas2 != 0.0) {
	norm = 2.0 / f_phas2;
	for (it=1; it<=ds->nphase; it++)
	  PHAS2(it,lu) *= norm;
      }

    } /* end for lu<ntau */

    prep_double_scat_integr (ds->nphase, ds->ntau, nf, ds->mu_phase,
			     phas2, mu_eq, neg_phas, norm_phas);
  } /* end if (need_secondary_scattering) */

  /* end of BDE stuff */

  /*
   * Start loop over zenith angles
   */
  for (iu = 1; iu <= ds->numu; iu++) {
    if (UMU(iu) < 0.) {
      /*
       * Calculate zenith angles of incident and emerging directions
       */
      theta0 = acos(-ds->bc.umu0)/DEG;
      thetap = acos(UMU(iu))/DEG;
    }
    /*
     * Start loop over azimuth angles
     */
    for (jp = 1; jp <= ds->nphi; jp++) {
      /*
       * Calculate cosine of scattering angle, eq. STWL(4)
       */
      ctheta = -ds->bc.umu0*UMU(iu)+sqrt((1.-SQR(ds->bc.umu0))*(1.-SQR(UMU(iu))))*cos(PHIRAD(jp));
      /*
       * Initialize phase function
       */
      for (lc = 1; lc <= ncut; lc++) {
        PHASM(lc) = 1.;
      }

      /* BDE ** Interpolate original phase function */
      /* BDE ** to actual phase function            */

      /* !!! +1: locate starts counting from 0! */
      it = locate_disort ( ds->mu_phase, ds->nphase, ctheta ) + 1;

      for (lc=1; lc<=ncut; lc++)
	PHASA(lc) = DSPHASE(it,lc)
	  + ( ctheta - ds->MUP(it) ) /
	  ( ds->MUP(it+1) - ds->MUP(it) ) *
	  ( DSPHASE(it+1,lc) - DSPHASE(it,lc) );
      /*
       * Initialize Legendre poly. recurrence
       */
      plm1 = 1.;
      plm2 = 0.;
      for (k = 1; k <= ds->nstr-1; k++) {
        /*
         * Calculate Legendre polynomial of P-sub-l by upward recurrence
         */
        pl   = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
        plm2 = plm1;
        plm1 = pl;

        /*
         * Calculate delta-M transformed phase function
         */
	for (lc=1; lc <= ncut; lc++) {
	  PHASM(lc) += (double)(2*k+1)*pl*(PMOM(k,lc)-FLYR(lc))/(1.-FLYR(lc));
	}
      }
      /*
       * Apply TMS method, eq. STWL(68)
       */
      for (lc = 1; lc <= ncut; lc++) {
        PHAST(lc) = PHASA(lc)/(1.-FLYR(lc)*SSALB(lc));
      }
      for (lu = 1; lu <= ds->ntau; lu++) {
        if (!lyrcut || LAYRU(lu) < ncut) {
          ussndm        = c_single_scat(dither,LAYRU(lu),ncut,phast,ds->ssalb,taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          ussp          = c_single_scat(dither,LAYRU(lu),ncut,phasm,oprim,    taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          UU(iu,lu,jp) += ussndm-ussp;
        }
      }
      if (UMU(iu) < 0. && fabs(theta0-thetap) <= dtheta) {
        /*
         * Emerging direction is in the aureole (theta0 +/- dtheta).
         * Apply IMS method for correction of secondary scattering below top level.
         */
        ltau = 1;
        if (UTAU(1) <= dither) {
          ltau = 2;
        }
        for (lu = ltau; lu <= ds->ntau; lu++) {
          if(!lyrcut || LAYRU(lu) < ncut) {
            duims = c_new_secondary_scat(ds,iu,lu,it,ctheta,flyr,
					 LAYRU(lu),tauc,
					 nf,
					 phas2, mu_eq, neg_phas,
					 NORM_PHAS(lu));
	    UU(iu,lu,jp) -= duims;
          }
        }
      }
    } /* end loop over azimuth angles */
  } /* end loop over zenith angles */

  free(mu_eq); free(norm_phas); free(neg_phas);
  free(phas2); free(phasr);

  return;
}

/*============================= end of c_new_intensity_correction() =====*/
