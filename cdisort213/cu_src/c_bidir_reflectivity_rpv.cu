// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_bidir_reflectivity_rpv() ==============*/

/*
  Computes the Rahman, Pinty, Verstraete BRDF.  The incident
  and outgoing cosine zenith angles are MU1 and MU2, respectively,
  and the relative azimuthal angle is PHI.  In this case the incident
  direction is where the radiation is coming from, so MU1>0 and
  the hot spot is MU2=MU1 and PHI=180 (the azimuth convention is
  different from the original Frank Evans code).
  The reference is:
  Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere
  Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable
  With NOAA Advanced Very High Resolution Radiometer Data,
  J. Geophys. Res., 98, 20791-20801.

  Translated from fortran to C by Robert Buras; original name RPV_REFLECTION

  INPUT:

    rho0   :  BRDF rpv: rho0
    k      :  BRDF rpv: k
    theta  :  BRDF rpv: theta
    sigma  :  BRDF rpv snow: sigma
    t1     :  BRDF rpv snow: t1
    t2     :  BRDF rpv snow: t2
    scale  :  BRDF rpv: scale
    mu1    :  Cosine of angle of reflection (positive)
    mu2    :  Cosine of angle of incidence (positive)
    phi    :  Difference of azimuth angles of incidence and reflection
                 (radians)
    badmu  :  minimally allowed value for mu1 and mu2

  LOCAL VARIABLES:

    ans    :  Return value

   Called by- c_bidir_reflectivity
-------------------------------------------------------------------------*/

double c_bidir_reflectivity_rpv ( rpv_brdf_spec *brdf,
                                  double         mu1,
				  double         mu2,
				  double         phi,
				  double         badmu )
{
  double
    m, f, h, cosphi, sin1, sin2, cosg, tan1, tan2, capg,
    hspot, t, g;
  double ans;

  /* This function needs more checking; some constraints are
     required to avoid albedos larger than 1; in particular,
     the BDREF is limited to 5 times the hotspot value to
     avoid extremely large values at low polar angles */


  /* Azimuth convention different from Frank Evans:
     Here PHI=0 means the backward direction while
     while in DISORT PHI=0 means forward. */
  phi = M_PI - phi;

  /* Don't allow mu's smaller than BADMU because
     the albedo is larger than 1 for those */
  if ( badmu > 0.0 ) {
    if ( mu1 < badmu )
      mu1 = badmu;
    if ( mu2 < badmu )
      mu2 = badmu;
  }

  /* Hot spot */
  hspot = brdf->rho0 * ( pow ( 2.0 * mu1 * mu1 * mu1 , brdf->k - 1.0 ) *
		   ( 1.0 - brdf->theta ) / ( 1.0 + brdf->theta ) / ( 1.0 + brdf->theta )
		   *  ( 2.0 - brdf->rho0 )
		   + brdf->sigma / mu1 ) * ( brdf->t1 * exp ( M_PI * brdf->t2 ) + 1.0 );

  /* Hot spot region */
  /* is this bug??? phi <= 1e-4 would be more sensible ... RPB */
  if (phi == 1e-4 && mu1 == mu2)
    return hspot * brdf->scale;

  m = pow ( mu1 * mu2 * ( mu1 + mu2 ) , brdf->k - 1.0 );
  cosphi = cos(phi);
  sin1 = sqrt ( 1.0 - mu1 * mu1 );
  sin2 = sqrt ( 1.0 - mu2 * mu2 );
  cosg = mu1 * mu2 + sin1 * sin2 * cosphi;
  g = acos ( cosg );
  f = ( 1.0 - brdf->theta * brdf->theta ) /
    pow ( 1.0 + 2.0 * brdf->theta * cosg + brdf->theta * brdf->theta , 1.5);

  tan1 = sin1 / mu1;
  tan2 = sin2 / mu2;
  capg = sqrt( tan1 * tan1 + tan2 * tan2 - 2.0 * tan1 * tan2 * cosphi );
  h = 1.0 + ( 1.0 - brdf->rho0 ) / ( 1.0 + capg );
  t = 1.0 + brdf->t1 * exp ( brdf->t2 * ( M_PI - g ) );

  ans = brdf->rho0 * ( m * f * h + brdf->sigma / mu1 ) * t * brdf->scale;

 if (ans < 0.0)
   ans = 0.0;

 return ans;
}

/*============================= end of c_bidir_reflectivity_rpv() =======*/
