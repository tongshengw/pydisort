// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_albtrans_intensity() ===================*/

/*
   Computes azimuthally-averaged intensity at top and bottom of medium
   (related to albedo and transmission of medium by reciprocity principles;
   see Ref S2).  User polar angles are used as incident beam angles.
   (This is a very specializedversion of user_intensities)

   ** NOTE **  User input values of UMU (assumed positive) are temporarily in
               upper locations of  UMU  and corresponding negatives are in
               lower locations (this makes GU come out right); the contents
               of the temporary UMU array are:
                   -UMU(ds->numu),..., -UMU(1), UMU(1),..., UMU(ds->numu)

   I N P U T    V A R I A B L E S:

       ds     :  Disort state variables
       gu     :  Eigenvectors interpolated to user polar angles (i.e., g in eq. SC(1), STWL(31ab))
       kk     :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       ll     :  Constants of integration in eq. SC(1), obtained by solving scaled version of eq. SC(5);
                 exponential term of eq. SC(12) not included
       nn     :  Order of double-Gauss quadrature (NSTR/2)
       taucpr :  Cumulative optical depth (delta-M-scaled)

   O U T P U T    V A R I A B L E:

       out->u0u : Diffuse azimuthally-averaged intensity at top and bottom of medium (directly transmitted component,
                  corresponding to bndint in user_intensities, is omitted).

   I N T E R N A L    V A R I A B L E S:

       dtau   :  Optical depth of a computational layer
       palint :  Non-boundary-forced intensity component
       utaupr :  Optical depths of user output levels (delta-M scaled)
       wk     :  Scratch vector for saving 'EXP' evaluations
       All the exponential factors (i.e., exp1, expn,... etc.)
       come from the substitution of constants of integration in
       eq. SC(12) into eqs. S1(8-9).  All have negative arguments.

   Called by- c_albtrans
 -------------------------------------------------------------------*/

void c_albtrans_intensity(disort_state *ds,
			  disort_output *out,
                          double       *gu,
                          double       *kk,
                          double       *ll,
                          int           nn,
                          double       *taucpr,
                          double       *wk)
{
  register int
    iq,iu,iumax,iumin,lc,lu;
  double
    denom,dtau,exp1,exp2,expn,mu,palint,sgn,utaupr[2];

  UTAUPR(1) = 0.;
  UTAUPR(2) = TAUCPR(ds->nlyr);

  for (lu = 1; lu <= 2; lu++) {
    if (lu == 1) {
      iumin = ds->numu/2+1;
      iumax = ds->numu;
      sgn   = 1.;
    }
    else {
      iumin = 1;
      iumax = ds->numu/2;
      sgn   = -1.;
    }

    /*
     * Loop over polar angles at which albedos/transmissivities desired
     * ( upward angles at top boundary, downward angles at bottom )
     */
    for (iu = iumin; iu <= iumax; iu++) {
      mu = UMU(iu);
      /*
       * Integrate from top to bottom computational layer
       */
      palint = 0.;
      for (lc = 1; lc <= ds->nlyr; lc++) {
        dtau = TAUCPR(lc)-TAUCPR(lc-1);
        exp1 = exp((UTAUPR(lu)-TAUCPR(lc-1))/mu);
        exp2 = exp((UTAUPR(lu)-TAUCPR(lc  ))/mu);
        /*
         * KK is negative
         */
        for (iq = 1; iq <= nn; iq++) {
          WK(iq) = exp(KK(iq,lc)*dtau);
          denom  = 1.+mu*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            /*
             * L'Hospital limit
             */
            expn = dtau/mu*exp2;
          }
          else {
            expn = (exp1*WK(iq)-exp2)*sgn/denom;
          }
          palint += GU(iu,iq,lc)*LL(iq,lc)*expn;
        }
        /*
         * KK is positive
         */
        for (iq = nn+1; iq <= ds->nstr; iq++) {
          denom = 1.+mu*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            expn = -dtau/mu*exp1;
          }
          else {
            expn = (exp1-exp2*WK(ds->nstr+1-iq))*sgn/denom;
          }
          palint += GU(iu,iq,lc)*LL(iq,lc)*expn;
        }
      }
      U0U(iu,lu) = palint;
    }
  }

  return;
}

/*============================= end of c_albtrans_intensity() ============*/
