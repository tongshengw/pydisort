// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_interp_coefficients_beam_source() =======*/

/*
     Find coefficients at user angle, necessary for later use in
     c_interp_source()
*/

/*

    I N P U T      V A R I A B L E S:

       cmu    :   Computational polar angles
       chtau  :   The optical depth in spherical geometry.
       delmo  :   Kronecker delta, delta-sub-m0
       fbeam  :   incident beam radiation at top
       gl     :   Phase function Legendre coefficients multiplied by (2l+1) and single-scatter albedo
       lc:    :   layer index
       mazim  :   order of azimuthal component
       nstr   :   number of streams
       numu   :   number of user angles
       taucpr :   delta-m-scaled optical depth
       xba    :   alfa in eq. KS(7)
       ylmu   :   Normalized associated Legendre polynomial at the user angles -umu-
       ylm0   :   Normalized associated Legendre polynomial at the beam angle

    O U T P U T     V A R I A B L E S:

       zb0u   :   x-sub-zero in KS(7) at user angles -umu-
       zb1u   :   x-sub-one in KS(7) at user angles -umu-
       zju    :  Solution vector Z-sub-zero after solving eq. SS(19), STWL(24b), at user angles -umu-

   Called by- c_disort

*/

void c_interp_coefficients_beam_source(disort_state   *ds,
				       double         *chtau,
				       double          delm0,
				       double          fbeam,
				       double         *gl,
				       int             lc,
				       int             mazim,
				       int             nstr,
				       int             numu,
				       double         *taucpr,
				       disort_triplet *zbu,
				       double         *xba,
				       double         *zju,
				       double         *ylm0,
				       double         *ylmu)
{
  register int
    iu,k;
  double
    deltat,sum,q0a,q2a,q0,q2;

  /*     Calculate x-sub-zero in STWJ(6d) */
  deltat = TAUCPR(lc) - TAUCPR(lc-1);

  q0a = exp(-CHTAU(lc-1));
  q2a = exp(-CHTAU(lc));

  for (iu = 1; iu <= numu; iu++) {
    sum = 0.0;
    for (k = mazim; k <= nstr-1; k++) {
      sum = sum + GL(k,lc)*YLMU(k,iu)*YLM0(k);
    }
    ZJU(iu) = (2.0-delm0)*fbeam*sum/(4.0*M_PI);
  }

  for (iu = 1; iu <= numu; iu++) {

    q0 = q0a*ZJU(iu);
    q2 = q2a*ZJU(iu);

    /*     x-sub-zero and x-sub-one in Eqs. KS(48-49)   */

    ZB1U(iu,lc)=(1./deltat)*(q2*exp(XBA(lc)*TAUCPR(lc))
			     -q0*exp(XBA(lc)*TAUCPR(lc-1)));
    ZB0U(iu,lc) = q0*exp(XBA(lc)*TAUCPR(lc-1))-ZB1U(iu,lc)*TAUCPR(lc-1);
  }

  return;

}
/*============================= end of c_interp_coefficients_beam_source() =*/
