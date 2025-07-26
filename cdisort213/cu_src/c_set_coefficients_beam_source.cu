// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_set_coefficients_beam_source() ========*/

/*
       Set coefficients in ks(7) for beam source

    I N P U T      V A R I A B L E S:

       cmu    :   Computational polar angles
       ch     :   The Chapman-factor to correct for pseudo-spherical geometry in the direct beam.
       chtau  :   The optical depth in spherical geometry.
       delmo  :   Kronecker delta, delta-sub-m0
       fbeam  :   incident beam radiation at top
       gl     :   Phase function Legendre coefficients multiplied by (2l+1) and single-scatter albedo
       lc:    :   layer index
       mazim  :   order of azimuthal component
       nstr   :   number of streams
       taucpr :   delta-m-scaled optical depth
       ylmc   :   Normalized associated Legendre polynomial at the quadrature angles -cmu-
       ylm0   :   Normalized associated Legendre polynomial at the beam angle

    O U T P U T     V A R I A B L E S:

       xba    :   alfa in eq. KS(7)
       xb0    :   x-sub-zero in KS(7)
       xb1    :   x-sub-one in KS(7)
       zj     :  Solution vector Z-sub-zero after solving eq. SS(19), STWL(24b)

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_set_coefficients_beam_source(disort_state *ds,
				    double       *ch,
				    double       *chtau,
				    double       *cmu,
				    double        delm0,
				    double        fbeam,
				    double       *gl,
				    int           lc,
				    int           mazim,
				    int           nstr,
				    double       *taucpr,
				    double       *xba,
				    disort_pair  *xb,
				    double       *ylm0,
				    double       *ylmc,
				    double       *zj)
{

  register int
    iq,k;
  double
    deltat,sum,q0a,q2a,q0,q2;
  static double
    big;

  big    = sqrt(DBL_MAX)/1.e+10;

  /*     Calculate x-sub-zero in STWJ(6d)   */

  for (iq = 1; iq <= nstr; iq++) {
    sum = 0;
    for (k = mazim; k <= nstr-1; k++) {
      sum += GL(k,lc)*YLMC(k,iq)*YLM0(k);
    }
    ZJ(iq) = (2.-delm0)*fbeam*sum/(4.*M_PI);
  }

  q0a = exp( -CHTAU(lc-1) );
  q2a = exp( -CHTAU(lc) );

  /*     Calculate alfa coefficient  */

  deltat = TAUCPR(lc) - TAUCPR(lc-1);

  XBA(lc) = 1./CH(lc);

  if ( fabs(XBA(lc)) > big  &&  TAUCPR(lc) > 1.)  XBA(lc) = 0.0;

  if( fabs(XBA(lc)*TAUCPR(lc)) > log(big))	  XBA(lc) = 0.0;

  /*     Dither alfa if it is close to one of the quadrature angles */

  if (  fabs(XBA(lc)) > 0.00001 ) {
    for (iq = 1; iq <= nstr/2; iq++) {
      if (fabs((fabs(XBA(lc))-1.0/CMU(iq))/XBA(lc) ) < 0.05 ) XBA(lc) = XBA(lc) * 1.001;
    }
  }

  for (iq = 1; iq <= nstr; iq++) {

    q0 = q0a * ZJ(iq);
    q2 = q2a * ZJ(iq);

    /*     x-sub-zero and x-sub-one in Eqs. KS(48-49)   */

    XB1(iq,lc) = (1.0/deltat)*(q2*exp(XBA(lc)*TAUCPR(lc)) - q0*exp(XBA(lc)*TAUCPR(lc-1)));
    XB0(iq,lc) = q0 * exp(XBA(lc)*TAUCPR(lc-1)) - XB1(iq,lc)*TAUCPR(lc-1);

  }
  return;
}
/*============================= end of c_set_coefficients_beam_source() ====*/
