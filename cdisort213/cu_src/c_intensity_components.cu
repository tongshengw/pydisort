// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_intensity_components() ================*/

/*
    Calculates the Fourier intensity components at the quadrature
    angles for azimuthal expansion terms (mazim) in eq. SD(2),STWL(6)

    I N P U T    V A R I A B L E S:

       ds      :  Disort state variables
       kk      :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       gc      :  Eigenvectors at polar quadrature angles in eq. SC(1)
       ll      :  Constants of integration in eq. SC(1), obtained by solving scaled version of eq. SC(5);
                  exponential term of eq. SC(12) not included
       lyrcut  :  Logical flag for truncation of computational layer
       mazim   :  Order of azimuthal component
       ncut    :  Number of computational layer where absorption optical depth exceeds ABSCUT
       nn      :  Order of double-Gauss quadrature (NSTR/2)
       taucpr  :  Cumulative optical depth (delta-M-scaled)
       utaupr  :  Optical depths of user output levels in delta-M coordinates;  equal to UTAU if no delta-M
       zz      :  Beam source vectors in eq. SS(19), STWL(24b)
       plk     :  Thermal source vectors z0,z1 by solving eq. SS(16), Y-sub-zero, Y-sub-one in STWL(26ab);
                  plk[].zero, plk[].one (see cdisort.h)

    O U T P U T   V A R I A B L E S:

       uum     :  Fourier components of the intensity in eq. SD(12) (at polar quadrature angles)

    I N T E R N A L   V A R I A B L E S:

       fact    :  exp(-utaupr/umu0)
       zint    :  intensity of m=0 case, in eq. SC(1)

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_intensity_components(disort_state *ds,
                            double       *gc,
                            double       *kk,
                            int          *layru,
                            double       *ll,
                            int           lyrcut,
                            int           mazim,
                            int           ncut,
                            int           nn,
                            double       *taucpr,
                            double       *utaupr,
                            double       *zz,
                            disort_pair  *plk,
                            double       *uum)
{
  register int
    iq,jq,lu,lyu;
  register double
    zint;

  /*
   * Loop over user levels
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    lyu = LAYRU(lu);
    if (lyrcut && lyu > ncut) {
      continue;
    }
    for (iq = 1; iq <= ds->nstr; iq++) {
      zint = 0.;
      for (jq = 1; jq <= nn; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu  )));
      }
      for (jq = nn+1; jq <=ds->nstr; jq++) {
        zint += GC(iq,jq,lyu)*LL(jq,lyu)*exp(-KK(jq,lyu)*(UTAUPR(lu)-TAUCPR(lyu-1)));
      }
      UUM(iq,lu) = zint;
      if (ds->bc.fbeam > 0.) {
        UUM(iq,lu) = zint+ZZ(iq,lyu)*exp(-UTAUPR(lu)/ds->bc.umu0);
      }
      if (ds->flag.planck && mazim == 0) {
        UUM(iq,lu) += ZPLK0(iq,lyu)+ZPLK1(iq,lyu)*UTAUPR(lu);
      }
    }
  }

  return;
}

/*============================= end of c_intensity_components() =========*/
