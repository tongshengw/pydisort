// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_single_scat() =========================*/

/*
        Calculates single-scattered intensity from eqs. STWL (65b,d,e)

                I N P U T   V A R I A B L E S

        dither   small multiple of machine precision
        layru    index of utau in multi-layered system
        nlyr     number of sublayers
        phase    phase functions of sublayers
        omega    single scattering albedos of sublayers
        tau      optical thicknesses of sublayers
        umu      cosine of emergent angle
        umu0     cosine of incident zenith angle
        utau     user defined optical depth for output intensity
        fbeam   incident beam radiation at top


   Called by- c_intensity_correction
 -------------------------------------------------------------------*/

double c_single_scat(double   dither,
                     int      layru,
                     int      nlyr,
                     double  *phase,
                     double  *omega,
                     double  *tau,
                     double   umu,
                     double   umu0,
                     double   utau,
                     double   fbeam)
{
  register int
    lyr;
  double
    ans,exp0,exp1;

  ans  = 0.;
  exp0 = exp(-utau/umu0);

  if (fabs(umu+umu0) <= dither) {
    /*
     * Calculate downward intensity when umu=umu0, eq. STWL (65e)
     */
    for (lyr = 1; lyr <= layru-1; lyr++) {
      ans += OMEGA(lyr)*PHASE(lyr)*(TAU(lyr)-TAU(lyr-1));
    }
    ans = fbeam/(4.*M_PI*umu0)*exp0*(ans+OMEGA(layru)*PHASE(layru)*(utau-TAU(layru-1)));
    return ans;
  }

  if (umu > 0.) {
    /*
     * Upward intensity, eq. STWL (65b)
     */
    for (lyr = layru; lyr <= nlyr; lyr++) {
      exp1  = exp(-((TAU(lyr)-utau)/umu+TAU(lyr)/umu0));
      ans  += OMEGA(lyr)*PHASE(lyr)*(exp0-exp1);
      exp0  = exp1;
    }
  }
  else {
    /*
     * Downward intensity, eq. STWL (65d)
     */
    for (lyr = layru; lyr >= 1; lyr--) {
      exp1  = exp(-((TAU(lyr-1)-utau)/umu+TAU(lyr-1)/umu0));
      ans  += OMEGA(lyr)*PHASE(lyr)*(exp0-exp1);
      exp0  = exp1;
    }
  }
  ans *= fbeam/(4.*M_PI*(1.+umu/umu0));

  return ans;
}

/*============================= end of c_single_scat() ==================*/
