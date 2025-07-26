// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_solve1() ===============================*/

/*
     Construct right-hand side vector -b- for isotropic incidence
     (only) on either top or bottom boundary and solve system
     of equations obtained from the boundary conditions and the
     continuity-of-intensity-at-layer-interface equations

     I N P U T      V A R I A B L E S:

       ds       :  Disort state variables
       cband    :  Left-hand side matrix of banded linear system
                   eq. SC(5), scaled by eq. SC(12); assumed already
                   in LU-decomposed form, ready for LINPACK solver
       ihom     :  Direction-of-illumination flag (TOP_ILLUM, top; BOT_ILLUM, bottom)
       ipvt     :
       ncol     :  Number of columns in CBAND
       ncut     :
       nn       :  Order of double-Gauss quadrature (NSTR/2)

    O U T P U T     V A R I A B L E S:

       b        :  Right-hand side vector of eq. SC(5) going into
                   sgbsl; returns as solution vector of eq.
                   SC(12), constants of integration without
                   exponential term
       ll       :  permanent storage for -b-, but re-ordered


    I N T E R N A L    V A R I A B L E S:

       ipvt     :  INTEGER vector of pivot indices
       ncd      :  Number of diagonals below or above main diagonal

   Called by- c_albtrans
   Calls- c_sgbsl
 +-------------------------------------------------------------------+
*/

void c_solve1(disort_state *ds,
              double       *cband,
              int           ihom,
              int          *ipvt,
              int           ncol,
              int           ncut,
              int           nn,
              double       *b,
              double       *ll)
{
  register int
    i,ipnt,iq,lc,ncd;

  memset(b,0,ds->nstr*ds->nlyr*sizeof(double));

  if (ihom == TOP_ILLUM) {
    /*
     * Because there are no beam or emission sources, remainder of B array is zero
     */
    for (i = 1; i <= nn; i++) {
      B(i)         = ds->bc.fisot;
      B(ncol-nn+i) = 0.;
    }
  }
  else if (ihom == BOT_ILLUM) {
    for (i = 1; i <= nn; i++) {
      B(i)         = 0.;
      B(ncol-nn+i) = ds->bc.fisot;
    }
  }
  else {
    c_errmsg("solve1---unrecognized ihom",DS_ERROR);
  }

  ncd = 3*nn-1;
  c_sgbsl(cband,(9*(ds->nstr/2)-2),ncol,ncd,ncd,ipvt,b,0);
  for (lc = 1; lc <= ncut; lc++) {
    ipnt = lc*ds->nstr-nn;
    for (iq = 1; iq <= nn; iq++) {
      LL(nn-iq+1,lc) = B(ipnt-iq+1);
      LL(nn+iq,  lc) = B(ipnt+iq  );
    }
  }

  return;
}

/*============================= end of c_solve1() ========================*/
