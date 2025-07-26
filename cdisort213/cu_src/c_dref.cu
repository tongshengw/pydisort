// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_dref() ================================*/

/*
  Flux albedo for given angle of incidence, given a bidirectional reflectivity.

  INPUTS
    wvnmlo    :  Lower wavenumber (inv-cm) of spectral interval
    wvnmhi    :  Upper wavenumber (inv-cm) of spectral interval
    mu        :  Cosine of incidence angle
    brdf_type :  BRDF type
    brdf      :  pointer to disort_brdf structure
    callnum   :  number of surface calls

  INTERNAL VARIABLES

    gmu    : The NMUG angle cosine quadrature points on (0,1)
             NMUG is set in cdisort.h
    gwt    : The NMUG angle cosine quadrature weights on (0,1)

   Called by- c_check_inputs
   Calls- c_gaussian_quadrature, c_errmsg, c_bidir_reflectivity
 --------------------------------------------------------------------*/

double c_dref(double       wvnmlo,
              double       wvnmhi,
              double       mu,
	      int          brdf_type,
	      disort_brdf *brdf,
	      int          callnum )
{
  static int
    pass1 = TRUE;
  register int
    jg,k;
  double
    ans,sum;
  static double
    gmu[NMUG],gwt[NMUG];

  if (pass1) {
    pass1 = FALSE;
    c_gaussian_quadrature(NMUG/2,gmu,gwt);
    for (k = 1; k <= NMUG/2; k++) {
      GMU(k+NMUG/2) = -GMU(k);
      GWT(k+NMUG/2) =  GWT(k);
    }
  }

  if (fabs(mu) > 1.) {
    c_errmsg("dref--input argument error(s)",DS_ERROR);
  }

  ans = 0.;
  /*
   * Loop over azimuth angle difference
   */
  for (jg = 1; jg <= NMUG; jg++) {
    /*
     * Loop over angle of reflection
     */
    sum = 0.;
    for (k = 1; k <= NMUG/2; k++) {
      sum += GWT(k) * GMU(k) *
	c_bidir_reflectivity ( wvnmlo, wvnmhi, GMU(k), mu, M_PI*GMU(jg), brdf_type, brdf, callnum );
    }
    ans += GWT(jg)*sum;
  }
  if (ans < 0. || ans > 1.) {
    c_errmsg("DREF--albedo value not in [0,1]",DS_WARNING);
  }

  return ans;
}

/*============================= end of c_dref() =========================*/
