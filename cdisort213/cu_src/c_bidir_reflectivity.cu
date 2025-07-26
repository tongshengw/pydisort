// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_bidir_reflectivity() ==================*/

/*
  Supplies surface bi-directional reflectivity.

  NOTE 1: Bidirectional reflectivity in DISORT is defined by eq. 39 in STWL.
  NOTE 2: Both MU and MU0 (cosines of reflection and incidence angles) are positive.

  Translated from fortran to C by Robert Buras; original name BDREF

  INPUT:

    wvnmlo    : Lower wavenumber (inv cm) of spectral interval
    wvnmhi    : Upper wavenumber (inv cm) of spectral interval
    mu        : Cosine of angle of reflection (positive)
    mup       : Cosine of angle of incidence (positive)
    dphi      : Difference of azimuth angles of incidence and reflection
                (radians)
    brdf_type : BRDF type
    brdf      : BRDF input
    callnum   : number of surface calls

  LOCAL VARIABLES:

    ans       :  Return variable
    badmu     :  minimally allowed value for mu1 and mu2
    flxalb    :
    irmu      :
    rmu       :
    swvnmlo   : value of wvnmlo from last call of this routine
    swvnmhi   : value of wvnmhi from last call of this routine
    srho0     : value of rho0   from last call of this routine
    sk        : value of k      from last call of this routine
    stheta    : value of theta  from last call of this routine
    ssigma    : value of sigma  from last call of this routine
    st1       : value of t1     from last call of this routine
    st2       : value of t2     from last call of this routine
    sscale    : value of scale  from last call of this routine
    siso      : value of iso    from last call of this routine
    svol      : value of vol    from last call of this routine
    sgeo      : value of geo    from last call of this routine
    su10      : value of u10    from last call of this routine
    spcl      : value of pcl    from last call of this routine
    ssal      : value of sal    from last call of this routine

   Called by- c_dref, c_surface_bidir
   Calls- c_dref, c_bidir_reflectivity_hapke,
          c_bidir_reflectivity_rpv, ocean_brdf, ambrals_brdf
-------------------------------------------------------------------------*/

double c_bidir_reflectivity ( double       wvnmlo,
			      double       wvnmhi,
			      double       mu,
			      double       mup,
			      double       dphi,
			      int          brdf_type,
			      disort_brdf *brdf,
			      int          callnum )
{
  int
    irmu;

  double
    ans, rmu, flxalb;

  static double
    badmu, swvnmlo, swvnmhi, srho0, sk,
    stheta, ssigma, st1, st2, sscale;

#if HAVE_BRDF
    static double
    siso, svol, sgeo;
#endif

  ans = 0.0;

  switch (brdf_type) {
  case BRDF_HAPKE:

    ans = c_bidir_reflectivity_hapke ( wvnmlo, wvnmhi, mu, mup, dphi );

    break;
  case BRDF_RPV:
    if ( swvnmlo != wvnmlo      ||
	 swvnmhi != wvnmhi      ||
	 srho0   != brdf->rpv->rho0   ||
	 sk      != brdf->rpv->k      ||
	 stheta  != brdf->rpv->theta  ||
	 ssigma  != brdf->rpv->sigma  ||
	 st1     != brdf->rpv->t1     ||
	 st2     != brdf->rpv->t2     ||
	 sscale  != brdf->rpv->scale ) {

      swvnmlo = wvnmlo;
      swvnmhi = wvnmhi;
      srho0   = brdf->rpv->rho0;
      sk      = brdf->rpv->k;
      stheta  = brdf->rpv->theta;
      ssigma  = brdf->rpv->sigma;
      st1     = brdf->rpv->t1;
      st2     = brdf->rpv->t2;
      sscale  = brdf->rpv->scale;

      badmu = 0.0;

      for (irmu=100; irmu>=0; irmu--) {

	rmu = ((double)irmu) * 0.01;

	flxalb = c_dref( wvnmlo, wvnmhi, rmu, brdf_type, brdf, callnum );

	if ( flxalb < 0.0 || flxalb > 1.0 ) {
	  badmu = rmu + 0.01;
	  if (badmu > 1.0)
	    badmu = 1.0;
	  printf("Using %f as limiting mu in RPV \n",badmu);
	  break;
	}
      }
    }

    ans = c_bidir_reflectivity_rpv ( brdf->rpv, mup, mu, dphi, badmu );

    break;
  case BRDF_CAM:

#if HAVE_BRDF
    /* call C tree saving function */
    /*
     * NOTE: Should group brdf->cam input arguments into the single pointer brdf->cam,
     *       in the same manner as brdf->rpv for c_bidir_reflectivity_rpv().
     */
    ans = ocean_brdf ( wvnmlo, wvnmhi, mu, mup, dphi,
		       brdf->cam->u10, brdf->cam->pcl, brdf->cam->xsal, callnum);

    /* remove BRDFs smaller than 0 */
    if (ans < 0.0)
      ans = 0.0;

    /* check for NaN */
    if ( ans != ans ) {
      printf("NaN returned from ocean_brdf: %e %e %e %e %e %e %e %e\n",	      wvnmlo, wvnmhi, mu, mup, dphi, brdf->cam->u10, brdf->cam->pcl, brdf->cam->xsal);
      ans = 1.0;
    }
#else
    c_errmsg("Error, ocean_brdf is not linked with your code!",DS_ERROR);
#endif
    break;
  case BRDF_AMB:

#if HAVE_BRDF
    /* mu = 0 or dmu = 0 cause problems */
    if ( siso != brdf->ambrals->iso ||
	 svol != brdf->ambrals->vol ||
	 sgeo != brdf->ambrals->geo ) {

      siso = brdf->ambrals->iso;
      svol = brdf->ambrals->vol;
      sgeo = brdf->ambrals->geo;

      badmu = 0.0;

      for (irmu=100; irmu>=0; irmu--) {

	rmu = ((double)irmu) * 0.01;

	flxalb = c_dref( wvnmlo, wvnmhi, rmu, brdf_type, brdf, callnum );

	if ( flxalb < 0.0 || flxalb > 1.0 ) {
	  badmu = rmu + 0.01;
	  if (badmu > 1.0)
	    badmu = 1.0;
	  printf("Using %f as limiting mu in AMBRALS \n",badmu);
	  break;
	}
      }
    }

    /* convert phi to degrees */
    /*    sdphi = dphi;
	  smup  = mup;
	  smu   = mu; probably no longer needed */

    dphi /= DEG;

    if ( badmu > 0.0 ) {
      if ( mu < badmu )
	mu = badmu;
      if ( mup < badmu )
	mup = badmu;
    }

    /*
     * NOTE: Should group brdf->ambrals input arguments into the single pointer brdf->ambrals,
     *       in the same manner as brdf->rpv for c_bidir_reflectivity_rpv().
     */
    ans = ambrals_brdf (brdf->ambrals->iso, brdf->ambrals->vol, brdf->ambrals->geo, mu, mup, dphi);

    /*    dphi = sdphi;
	  mup  = smup;
	  mu   = smu; probably no longer needed */

    /* check for NaN */
    if ( ans != ans ) {
      printf("NaN returned from ambrals_brdf: %e %e %e %e %e %e %e %e\n",	      wvnmlo, wvnmhi, mu, mup, dphi, brdf->ambrals->iso, brdf->ambrals->vol, brdf->ambrals->geo);
      ans = 1.0;
    }
#else
    c_errmsg("Error, ambrals_brdf is not linked with your code!",DS_ERROR);
#endif

    break;
  default:
    printf("bidir_reflectivity--surface BDRF model %d not known",	    brdf_type);
    c_errmsg("Exiting...",DS_ERROR);
  }

  return ans;
}

/*============================= end of c_bidir_reflectivity() ===========*/
