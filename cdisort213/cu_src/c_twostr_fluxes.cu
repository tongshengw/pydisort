// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_fluxes() ========================*/

/*
 Calculates the radiative fluxes, mean intensity, and flux derivative
 with respect to optical depth from the azimuthally-averaged intensity

 I n p u t     v a r i a b l e s:

   ds         :  'Disort' state variables
   ts         :  twostr_xyz structure variables (xp_0, yb_0d, zb_a...; see cdisort.h)
   ch         :  Chapman factor
   cmu        :  Abscissa for gauss quadrature over angle cosine
   kk         :  Eigenvalues
   layru      :  Layer numbers of user levels -utau-
   ll         :  Constants of integration in eqs. KST(42-43), obtaine by solving eqs. KST(38-41)
   lyrcut     :  Logical flag for truncation of comput. layer
   ncut       :  Number of computational layer where absorption optical depth exceeds -abscut-
   oprim      :  Delta-m scaled single scattering albedo
   rr         :  Eigenvectors at polar quadrature angles
   flag.spher :  TRUE turns on pseudo-spherical effects
   taucpr     :  Cumulative optical depth (delta-m-scaled)
   utaupr     :  Optical depths of user output levels in delta-m coordinates; equal to  -utau- if no delta-m

 O u t p u t     v a r i a b l e s:

   out      :  'Disort' output variables
   u0c      :  Azimuthally averaged intensities at polar quadrature angle cmu

 I n t e r n a l       v a r i a b l e s:

   dirint   :  direct intensity attenuated
   fdntot   :  total downward flux (direct + diffuse)
   fldir    :  fl[].zero, direct-beam flux (delta-m scaled)
   fldn     :  fl[].one, diffuse down-flux (delta-m scaled)
   fnet     :  net flux (total-down - diffuse-up)
   fact     :  EXP( - utaupr / ch ), where ch is the Chapman factor
   plsorc   :  Planck source function (thermal)
 ---------------------------------------------------------------------*/

void c_twostr_fluxes(disort_state  *ds,
                     twostr_xyz    *ts,
                     double        *ch,
                     double         cmu,
                     double        *kk,
                     int           *layru,
                     double        *ll,
                     int            lyrcut,
                     int            ncut,
                     double        *oprim,
                     double        *rr,
                     double        *taucpr,
                     double        *utaupr,
                     disort_output *out,
                     double        *u0c,
                     disort_pair   *fl)
{
  register int
    lu,lyu;
  double
    fdntot,fnet,plsorc,dirint;
  register double
    fact1,fact2;

  if (ds->flag.prnt[1]) {
    printf("\n\n                     <----------------------- Fluxes ----------------------->\n"                   "   optical  compu    downward    downward    downward       upward                    mean      Planck   d(net flux)\n"
                   "     depth  layer      direct     diffuse       total      diffuse         net   intensity      source   / d(op dep)\n");
  }

  memset(out->rad,0,ds->ntau*sizeof(disort_radiant));

  /*
   * Loop over user levels
   */
  if (ds->flag.planck) {
    for (lu = 1; lu <= ds->ntau; lu++) {
      lyu        = LAYRU(lu);
      fact1      = exp(-ZP_A(lyu)*UTAUPR(lu));
      U0C(1,lu) += fact1*(YP_0D(lyu)+YP_1D(lyu)*UTAUPR(lu));
      U0C(2,lu) += fact1*(YP_0U(lyu)+YP_1U(lyu)*UTAUPR(lu));
    }
  }
  for (lu = 1; lu <= ds->ntau; lu++) {
    lyu = LAYRU(lu);
    if (lyrcut && lyu > ncut) {
      /*
       * No radiation reaches this level
       */
      fdntot = 0.;
      fnet   = 0.;
      plsorc = 0.;
    }
    else {
      if (ds->bc.fbeam > 0.) {
        fact1      = exp(-ZB_A(lyu)*UTAUPR(lu));
        U0C(1,lu) += fact1*(YB_0D(lyu)+YB_1D(lyu)*UTAUPR(lu));
        U0C(2,lu) += fact1*(YB_0U(lyu)+YB_1U(lyu)*UTAUPR(lu));
        if (ds->bc.umu0 > 0. || ds->flag.spher) {
          fact1      = ds->bc.fbeam*exp(-UTAUPR(lu)/CH(lyu));
          dirint     = fact1;
          FLDIR(lu)  = fabs(ds->bc.umu0)*fact1;
          RFLDIR(lu) = fabs(ds->bc.umu0)*ds->bc.fbeam*exp(-UTAU(lu)/CH(lyu));
        }
        else {
          dirint     = 0.;
          FLDIR(lu)  = 0.;
          RFLDIR(lu) = 0.;
        }
      }
      else {
        dirint     = 0.;
        FLDIR(lu)  = 0.;
        RFLDIR(lu) = 0.;
      }
      fact1      = LL(1,lyu)*exp( KK(lyu)*(UTAUPR(lu)-TAUCPR(lyu  )));
      fact2      = LL(2,lyu)*exp(-KK(lyu)*(UTAUPR(lu)-TAUCPR(lyu-1)));
      U0C(1,lu) += fact2+RR(lyu)*fact1;
      U0C(2,lu) += fact1+RR(lyu)*fact2;
      /*
       * Calculate fluxes and mean intensities; downward and upward fluxes from eq. KST(9)
       */
      fact1     = 2.*M_PI*cmu;
      FLDN(lu)  = fact1*U0C(1,lu);
      FLUP(lu)  = fact1*U0C(2,lu);
      fdntot    = FLDN(lu)+FLDIR(lu);
      fnet      = fdntot-FLUP(lu);
      RFLDN(lu) = fdntot-RFLDIR(lu);
      /*
       * Mean intensity from eq. KST(10)
       */
      UAVG(lu) = U0C(1,lu)+U0C(2,lu);
      UAVG(lu) = (2.*M_PI*UAVG(lu)+dirint)/(4.*M_PI);

      /*
       * Flux divergence from eqs. KST(11-12)
       */
      plsorc   = 1./(1.-OPRIM(lyu))*exp(-ZP_A(lyu)*UTAUPR(lu))*(XP_0(lyu)+XP_1(lyu)*UTAUPR(lu));
      DFDT(lu) = (1.-SSALB(lyu))*4.*M_PI*(UAVG(lu)-plsorc);
    }
    if (ds->flag.prnt[1]) {
      printf("%10.4f%7d%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%14.3e\n",                     UTAU(lu),lyu,RFLDIR(lu),RFLDN(lu),fdntot,FLUP(lu),fnet,UAVG(lu),plsorc,DFDT(lu));
    }
  }

  return;
}

/*============================= end of c_twostr_fluxes() =================*/
