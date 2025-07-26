// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_set_matrix() ==========================*/

/*
    Calculate coefficient matrix for the set of equations obtained from the
    boundary conditions and the continuity-of-intensity-at-layer-interface equations.

    Store in the special banded-matrix format required by LINPACK routines


    I N P U T      V A R I A B L E S:

       ds       :  Disort state variables
       bdr      :  surface bidirectional reflectivity
       cmu,cwt  :  abscissae, weights for Gauss quadrature over angle cosine
       delm0    :  Kronecker delta, delta-sub-m0
       gc       :  Eigenvectors at polar quadrature angles, SC(1)
       kk       :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       lyrcut   :  Logical flag for truncation of computational layers
       nn       :  Number of streams in a hemisphere (NSTR/2)
       ncut     :  Total number of computational layers considered
       taucpr   :  Cumulative optical depth (delta-M-scaled)

   O U T P U T     V A R I A B L E S:

       cband    :  Left-hand side matrix of linear system eq. SC(5), scaled by eq. SC(12);
                   in banded form required by LINPACK solution routines
       ncol     :  Number of columns in cband


   I N T E R N A L    V A R I A B L E S:

       irow     :  Points to row in CBAND
       jcol     :  Points to position in layer block
       lda      :  Row dimension of CBAND
       ncd      :  Number of diagonals below or above main diagonal
       nshift   :  For positioning number of rows in band storage
       wk       :  Temporary storage for EXP evaluations


   BAND STORAGE

      LINPACK requires band matrices to be input in a special
      form where the elements of each diagonal are moved up or
      down (in their column) so that each diagonal becomes a row.
      (The column locations of diagonal elements are unchanged.)

      Example:  if the original matrix is

          11 12 13  0  0  0
          21 22 23 24  0  0
           0 32 33 34 35  0
           0  0 43 44 45 46
           0  0  0 54 55 56
           0  0  0  0 65 66

      then its LINPACK input form would be:

           *  *  *  +  +  +  , * = not used
           *  * 13 24 35 46  , + = used for pivoting
           * 12 23 34 45 56
          11 22 33 44 55 66
          21 32 43 54 65  *

      If A is a band matrix, the following program segment
      will convert it to the form (ABD) required by LINPACK
      band-matrix routines:

        n  = (column dimension of a, abd)
        ml = (band width below the diagonal)
        mu = (band width above the diagonal)
        m = ml+mu+1;
        for (j = 1; j <= n; j++) {
          i1 = IMAX(1,j-mu);
          i2 = IMIN(n,j+ml);
          for (i = i1; i <= i2; i++) {
            k = i-j+m;
            ABD(k,j) = A(i,j);
          }
        }

      This uses rows  ml+1 through  2*ml+mu+1  of ABD.
      The total number of rows needed in ABD is 2*ml+mu+1.
      In the example above, n = 6, ml = 1, mu = 2, and the
      row dimension of ABD must be >= 5.

   Called by- c_disort, c_albtrans
 -------------------------------------------------------------------*/

void c_set_matrix(disort_state *ds,
                  double       *bdr,
                  double       *cband,
                  double       *cmu,
                  double       *cwt,
                  double        delm0,
                  double       *dtaucpr,
                  double       *gc,
                  double       *kk,
                  int           lyrcut,
                  int          *ncol,
                  int           ncut,
                  double       *taucpr,
                  double       *wk)
{
  int
    mi     = ds->nstr/2,
    mi9m2  = 9*mi-2,
    nnlyri = ds->nstr*ds->nlyr,
    nn     = ds->nstr/2;
  register int
    iq,irow,jcol,jq,k,lc,lda,ncd,nncol,nshift;
  double
    expa,sum;

  memset(cband,0,mi9m2*nnlyri*sizeof(double));

  ncd    = 3*nn-1;
  lda    = 3*ncd+1;
  nshift = lda-2*ds->nstr+1;
  *ncol  = 0;

  /*
   * Use continuity conditions of eq. STWJ(17) to form coefficient matrix in STWJ(20);
   * employ scaling transformation STWJ(22)
   */
  for (lc = 1; lc <= ncut; lc++) {
    for (iq = 1; iq <= nn; iq++) {
      WK(iq) = exp(KK(iq,lc)*DTAUCPR(lc));
    }
    jcol = 0;
    for (iq = 1; iq <= nn; iq++) {
      *ncol += 1;
      irow   = nshift-jcol;
      for (jq = 1; jq <= ds->nstr; jq++) {
        CBAND(irow+ds->nstr,*ncol) =  GC(jq,iq,lc);
        CBAND(irow,         *ncol) = -GC(jq,iq,lc)*WK(iq);
        irow++;
      }
      jcol++;
    }

    for (iq = nn+1; iq <= ds->nstr; iq++) {
      *ncol += 1;
      irow = nshift-jcol;
      for (jq = 1; jq <= ds->nstr; jq++) {
        CBAND(irow+ds->nstr,*ncol) =  GC(jq,iq,lc)*WK(ds->nstr+1-iq);
        CBAND(irow,         *ncol) = -GC(jq,iq,lc);
        irow++;
      }
      jcol++;
    }
  }

  /*
   * Use top boundary condition of STWJ(20a) for first layer
   */
  jcol = 0;
  for (iq = 1; iq <= nn; iq++) {
    expa = exp(KK(iq,1)*TAUCPR(1));
    irow = nshift-jcol+nn;
    for (jq = nn; jq >= 1; jq--) {
      CBAND(irow,jcol+1) = GC(jq,iq,1)*expa;
      irow++;
    }
    jcol++;
  }

  for (iq = nn+1; iq <=ds->nstr; iq++) {
    irow = nshift-jcol+nn;
    for (jq = nn; jq >= 1; jq--) {
      CBAND(irow,jcol+1) = GC(jq,iq,1);
      irow++;
    }
    jcol++;
  }

  /*
   * Use bottom boundary condition of STWJ(20c) for last layer
   */
  nncol = *ncol-ds->nstr;
  jcol  = 0;
  for (iq = 1; iq <= nn; iq++) {
    nncol++;
    irow = nshift-jcol+ds->nstr;
    for (jq = nn+1; jq <= ds->nstr; jq++) {
      if (lyrcut || ( ds->flag.lamber && delm0 == 0. ) ) {
        /*
         * No azimuthal-dependent intensity if Lambert surface;
         * no intensity component if truncated bottom layer
         */
        CBAND(irow,nncol) = GC(jq,iq,ncut);
      }
      else {
        sum = 0.;
        for (k = 1; k <= nn; k++) {
          sum += CWT(k)*CMU(k)*BDR(jq-nn,k)*GC(nn+1-k,iq,ncut);
        }
        CBAND(irow,nncol) = GC(jq,iq,ncut)-(1.+delm0)*sum;
      }
      irow++;
    }
    jcol++;
  }

  for (iq = nn+1; iq <= ds->nstr; iq++) {
    nncol++;
    irow = nshift-jcol+ds->nstr;
    expa = WK(ds->nstr+1-iq);
    for (jq = nn+1; jq <= ds->nstr; jq++) {
      if (lyrcut || (ds->flag.lamber && delm0 == 0.)) {
        CBAND(irow,nncol) = GC(jq,iq,ncut)*expa;
      }
      else {
        sum = 0.;
        for (k = 1; k <= nn; k++) {
          sum += CWT(k)*CMU(k)*BDR(jq-nn,k)*GC(nn+1-k,iq,ncut);
        }
        CBAND(irow,nncol) = (GC(jq,iq,ncut)-(1.+delm0)*sum)*expa;
      }
      irow++;
    }
    jcol++;
  }

  return;
}

/*============================= end of c_set_matrix() ===================*/
