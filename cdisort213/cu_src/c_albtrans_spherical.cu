// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_albtrans_spherical() ===================*/

/*
    Calculates spherical albedo and transmissivity for the entire medium
    from the m=0 intensity components (this is a specialized version of fluxes)

    I N P U T    V A R I A B L E S:

       ds      :  Disort state variables
       cmu,cwt :  Abscissae, weights for Gaussian quadrature over angle cosine
       kk      :  Eigenvalues of coeff. matrix in eq. SS(7)
       gc      :  Eigenvectors at polar quadrature angles, SC(1)
       ll      :  Constants of integration in eq. SC(1), obtained by solving
                  scaled version of eq. SC(5); exponential term of eq. SC(12) not incl.
       nn      :  Order of double-Gauss quadrature (NSTR/2)

    O U T P U T   V A R I A B L E S:

       sflup   :  Up-flux at top (equivalent to spherical albedo due to
                  reciprocity).  For illumination from below it gives
                  spherical transmissivity

       sfldn   :  Down-flux at bottom (for single layer, equivalent to
                  spherical transmissivity due to reciprocity)

    I N T E R N A L   V A R I A B L E S:

       zint    :  Intensity of m=0 case, in eq. SC(1)

   Called by- c_albtrans
 --------------------------------------------------------------------*/

void c_albtrans_spherical(disort_state *ds,
                          double       *cmu,
                          double       *cwt,
                          double       *gc,
                          double       *kk,
                          double       *ll,
                          int           nn,
                          double       *taucpr,
                          double       *sflup,
                          double       *sfldn)
{
  register int
    iq,jq;
  double
    zint;

  *sflup = 0.;
  for (iq = nn+1; iq <= ds->nstr; iq++) {
    zint = 0.;
    for (jq = 1; jq <= nn; jq++) {
      zint += GC(iq,jq,1)*LL(jq,1)*exp(KK(jq,1)*TAUCPR(1));
    }
    for (jq = nn+1; jq <= ds->nstr; jq++) {
      zint += GC(iq,jq,1)*LL(jq,1);
    }
    *sflup += CWT(iq-nn)*CMU(iq-nn)*zint;
  }

  *sfldn = 0.;
  for (iq = 1; iq <= nn; iq++) {
    zint = 0.;
    for (jq = 1; jq <= nn; jq++) {
      zint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr);
    }
    for (jq = nn+1; jq <=ds->nstr; jq++) {
      zint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr)*exp(-KK(jq,ds->nlyr)*(TAUCPR(ds->nlyr)-TAUCPR(ds->nlyr-1)));
    }
    *sfldn += CWT(nn+1-iq)*CMU(nn+1-iq)*zint;
  }

  *sflup *= 2.;
  *sfldn *= 2.;

  return;
}

/*============================= end of c_albtrans_spherical() ============*/
