// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_upisot() ==============================*/

/*
    Finds the particular solution of thermal radiation of STWL(25)

    I N P U T     V A R I A B L E S:

       ds     :  Disort state variables
       cc     :  C-sub-ij in eq. SS(5), STWL(8b)
       cmu    :  Abscissae for Gauss quadrature over angle cosine
       oprim  :  Delta-M scaled single scattering albedo
       xr     :  Expansion coefficient b-sub-zero, b-sub-one of thermal source function, eq. STWL(24c)

    O U T P U T    V A R I A B L E S:

       zee    :  Solution vectors Z-sub-zero, Z-sub-one of eq. SS(16), STWL(26a,b)
       plk    :  Permanent storage for zee, but re-ordered

   I N T E R N A L    V A R I A B L E S:

       array  :  Coefficient matrix in left-hand side of eq. SS(16)
       ipvt   :  Integer vector of pivot indices required by LINPACK
       wk     :  Scratch array required by LINPACK

   Called by- c_disort
   Calls- c_sgeco, c_errmsg, c_sgesl
 -------------------------------------------------------------------*/

#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*ds->nstr]

void c_upisot(disort_state *ds,
              int           lc,
              double       *array,
              double       *cc,
              double       *cmu,
              int          *ipvt,
              int           nn,
              double       *oprim,
              double       *wk,
              disort_pair  *xr,
              disort_pair  *zee,
              disort_pair  *plk)
{
  register int
    iq,jq;
  double
    rcond;

  for (iq = 1; iq <= ds->nstr; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      ARRAY(iq,jq) = -CC(iq,jq);
    }
    ARRAY(iq,iq) += 1.;
    Z1(iq) = (1.-OPRIM(lc))*XR1(lc);
  }
  /*
   * Solve linear equations: same as in upbeam, except zj replaced by z1 and z0
   */
  rcond = 0.;
  c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);

  if (1.+rcond == 1.) {
    c_errmsg("upisot--sgeco says matrix near singular",DS_WARNING);
  }

  for (iq = 1; iq <= ds->nstr; iq++) {
    /* Need to use WK() as a buffer, since Z1 is part of a structure */
    WK(iq) = Z1(iq);
  }
  c_sgesl(array,ds->nstr,ds->nstr,ipvt,wk,0);
  for (iq = 1; iq <= ds->nstr; iq++) {
    Z1(iq) = WK(iq);
  }

  for (iq = 1; iq <= ds->nstr; iq++) {
    Z0(iq) = (1.-OPRIM(lc))*XR0(lc)+CMU(iq)*Z1(iq);
  }

  for (iq = 1; iq <= ds->nstr; iq++) {
    /* Need to use WK() as a buffer, since Z0 is part of a structure */
    WK(iq) = Z0(iq);
  }
  c_sgesl(array,ds->nstr,ds->nstr,ipvt,wk,0);
  for (iq = 1; iq <= ds->nstr; iq++) {
    Z0(iq) = WK(iq);
  }
  for (iq = 1; iq <= nn; iq++) {
    ZPLK0(nn+iq,  lc) = Z0(iq   );
    ZPLK1(nn+iq,  lc) = Z1(iq   );
    ZPLK0(nn-iq+1,lc) = Z0(iq+nn);
    ZPLK1(nn-iq+1,lc) = Z1(iq+nn);
  }

  return;
}

/*============================= end of c_upisot() =======================*/
