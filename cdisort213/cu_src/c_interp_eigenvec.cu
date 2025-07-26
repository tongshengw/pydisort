// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_interp_eigenvec() =====================*/

/*
   Interpolate eigenvectors to user angles; eq SD(8)

   Called by- c_disort, c_albtrans
 --------------------------------------------------------------------*/

void c_interp_eigenvec(disort_state *ds,
                       int           lc,
                       double       *cwt,
                       double       *evecc,
                       double       *gl,
                       double       *gu,
                       int           mazim,
                       int           nn,
                       double       *wk,
                       double       *ylmc,
                       double       *ylmu)
{
  register int
    iq,iu,jq,l;
  double
    sum;

  for (iq = 1; iq <= ds->nstr; iq++) {
    for (l = mazim; l <= ds->nstr-1; l++) {
      /*
       * Inner sum in SD(8) times all factors in outer sum but PLM(mu)
       */
      sum = 0.;
      for (jq = 1; jq <= ds->nstr; jq++) {
        sum += CWT(jq)*YLMC(l,jq)*EVECC(jq,iq);
      }
      WK(l+1) = .5*GL(l,lc)*sum;
    }
    /*
     * Finish outer sum in SD(8) and store eigenvectors
     */
    for (iu = 1; iu <= ds->numu; iu++) {
      sum = 0.;
      for (l = mazim; l <=ds->nstr-1; l++) {
        sum += WK(l+1)*YLMU(l,iu);
      }
      if (iq <= nn) {
        GU(iu,nn+iq,lc) = sum;
      }
      else {
        GU(iu,ds->nstr+1-iq,lc) = sum;
      }
    }
  }

  return;
}

/*============================= end of c_interp_eigenvec() ==============*/
