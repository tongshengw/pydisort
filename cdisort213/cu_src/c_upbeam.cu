// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_upbeam() ==============================*/

/*
   Finds the incident-beam particular solution of SS(18), STWL(24a)

   I N P U T    V A R I A B L E S:

       ds     :  Disort state variables
       cc     :  C-sub-ij in eq. SS(5)
       cmu    :  Abscissae for Gauss quadrature over angle cosine
       delm0  :  Kronecker delta, delta-sub-m0
       gl     :  Delta-M scaled Legendre coefficients of phase function
                 (including factors 2l+1 and single-scatter albedo)
       mazim  :  Order of azimuthal component
       ylm0   :  Normalized associated Legendre polynomial at the beam angle
       ylmc   :  Normalized associated Legendre polynomial at the quadrature angles

   O U T P U T    V A R I A B L E S:

       zj     :  Right-hand side vector X-sub-zero in SS(19),STWL(24b);
                 also the solution vector Z-sub-zero after solving that system
       zz     :  Permanent storage for zj, but re-ordered

   I N T E R N A L    V A R I A B L E S:

       array  :  Coefficient matrix in left-hand side of eq. SS(19), STWL(24b)
       ipvt   :  Integer vector of pivot indices required by LINPACK
       wk     :  Scratch array required by LINPACK

   Called by- c_disort
   Calls- c_sgeco, c_errmsg, c_sgesl
 -------------------------------------------------------------------*/

#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*ds->nstr]

void c_upbeam(disort_state *ds,
              int           lc,
              double       *array,
              double       *cc,
              double       *cmu,
              double        delm0,
              double       *gl,
              int          *ipvt,
              int           mazim,
              int           nn,
              double       *wk,
              double       *ylm0,
              double       *ylmc,
              double       *zj,
              double       *zz)
{
  register int
    iq,jq,k;
  double
    rcond,sum;

  for (iq = 1; iq <= ds->nstr; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      ARRAY(iq,jq) = -CC(iq,jq);
    }
    ARRAY(iq,iq) += 1.+CMU(iq)/ds->bc.umu0;
    sum = 0.;
    for (k = mazim; k <=ds->nstr-1; k++) {
      sum += GL(k,lc)*YLMC(k,iq)*YLM0(k);
    }
    ZJ(iq) = (2.-delm0)*ds->bc.fbeam*sum/(4.*M_PI);
  }
  /*
   * Find L-U (lower/upper triangular) decomposition of ARRAY and see if it is nearly singular
   * (NOTE:  ARRAY is altered)
   */
  rcond = 0.;
  c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);

  if (1.+rcond == 1.) {
    c_errmsg("upbeam--sgeco says matrix near singular",DS_WARNING);
  }

  /*
   * Solve linear system with coeff matrix ARRAY (assumed already L-U decomposed) and R.H. side(s) ZJ;
   * return solution(s) in ZJ
   */
  c_sgesl(array,ds->nstr,ds->nstr,ipvt,zj,0);
  for (iq = 1; iq <= nn; iq++) {
    ZZ(nn+iq,  lc) = ZJ(iq);
    ZZ(nn-iq+1,lc) = ZJ(iq+nn);
  }

  return;
}

/*============================= end of c_upbeam() =======================*/
