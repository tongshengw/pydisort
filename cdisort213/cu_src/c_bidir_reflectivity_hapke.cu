// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_bidir_reflectivity_hapke() ============*/

/*
 * Hapke's BRDF model (times Pi/Mu0):
 *   Hapke, B., Theory of reflectance and emittance spectroscopy, Cambridge University Press, 1993,
 * eq. 8.89 on page 233. Parameters are from Fig. 8.15 on page 231, except for w.

  INPUT:

    wvnmlo : Lower wavenumber (inv cm) of spectral interval
    wvnmhi : Upper wavenumber (inv cm) of spectral interval
    mu     : Cosine of angle of reflection (positive)
    mup    : Cosine of angle of incidence (positive)
    dphi   : Difference of azimuth angles of incidence and reflection
                (radians)

  LOCAL VARIABLES:

    iref   : bidirectional reflectance options; 1 - Hapke's BDR model
    b0     : empirical factor to account for the finite size of particles in Hapke's BDR model
    b      : term that accounts for the opposition effect (retroreflectance, hot spot) in Hapke's BDR model
    ctheta : cosine of phase angle in Hapke's BDR model
    gamma  : albedo factor in Hapke's BDR model
    h0     : H(mu0) in Hapke's BDR model
    h      : H(mu) in Hapke's BDR model
    hh     : angular width parameter of opposition effect in Hapke's BDR model
    p      : scattering phase function in Hapke's BDR model
    theta  : phase angle (radians); the angle between incidence and reflection directions in Hapke's BDR model
    w      : single scattering albedo in Hapke's BDR model

   Called by- c_bidir_reflectivity
-------------------------------------------------------------------------*/

double c_bidir_reflectivity_hapke ( double wvnmlo,
				    double wvnmhi,
				    double mu,
				    double mup,
				    double dphi )
{
  double
    b0,b,ctheta,Xgamm,
    h0,h,hh,p,thetah,w;

  ctheta = mu*mup+sqrt((1.-mu*mu)*(1.-mup*mup))*cos(dphi);
  thetah = acos(ctheta);
  p      = 1.+.5*ctheta;
  hh     =  .06;
  b0     = 1.;
  b      = b0*hh/(hh+tan(.5*thetah));
  w      = 0.6;
  Xgamm  = sqrt(1.-w);
  h0     = (1.+2.*mup)/(1.+2.*Xgamm*mup);
  h      = (1.+2.*mu )/(1.+2.*Xgamm*mu );

  return .25*w*((1.+b)*p+h0*h-1.0)/(mu+mup);
}

/*============================= end of c_bidir_reflectivity_hapke() =====*/
