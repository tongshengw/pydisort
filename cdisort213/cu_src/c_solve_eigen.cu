// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_solve_eigen() =========================*/

/*
   Solves eigenvalue/vector problem necessary to construct homogeneous
   part of discrete ordinate solution; STWJ(8b), STWL(23f)
   ** NOTE ** Eigenvalue problem is degenerate when single scattering
              albedo = 1;  present way of doing it seems numerically more
              stable than alternative methods that we tried

   I N P U T     V A R I A B L E S:

       ds     :  Disort state variables
       lc     :
       gl     :  Delta-M scaled Legendre coefficients of phase function
                 (including factors 2l+1 and single-scatter albedo)
       cmu    :  Computational polar angle cosines
       cwt    :  Weights for quadrature over polar angle cosine
       mazim  :  Order of azimuthal component
       nn     :  Half the total number of streams
       ylmc   :  Normalized associated Legendre polynomial
                 at the quadrature angles CMU


   O U T P U T    V A R I A B L E S:

       cc     :  C-sub-ij in eq. SS(5); needed in SS(15&18)
       eval   :  NN eigenvalues of eq. SS(12), STWL(23f) on return
                 from asymmetric_matrix but then square roots taken
       evecc  :  NN eigenvectors  (G+) - (G-)  on return
                 from asymmetric_matrix ( column j corresponds to EVAL(j) )
                 but then  (G+) + (G-)  is calculated from SS(10),
                 G+  and  G-  are separated, and  G+  is stacked on
                 top of  G-  to form NSTR eigenvectors of SS(7)
       gc     :  Permanent storage for all NSTR eigenvectors, but
                 in an order corresponding to KK
       kk     :  Permanent storage for all NSTR eigenvalues of SS(7),
                 but re-ordered with negative values first ( square
                 roots of EVAL taken and negatives added )


   I N T E R N A L   V A R I A B L E S:

       ab            :  Matrices AMB (alpha-beta), APB (alpha+beta) in reduced eigenvalue problem (see cdisort.h)
       array         :  Complete coefficient matrix of reduced eigenvalue
                        problem: (alpha+beta)*(alpha-beta)
       gpplgm        :  (g+) + (g-) (cf. eqs. SS(10-11))
       gpmigm        :  (g+) - (g-) (cf. eqs. SS(10-11))
       wk            :  Scratch array required by asymmetric_matrix

   Called by- c_disort, c_albtrans
   Calls- c_asymmetric_matrix, c_errmsg
 -------------------------------------------------------------------*/

/*
 * NOTE: Here the scratch array ARRAY(,) is half the size in each dimension compared to other subroutines
 */
#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*(ds->nstr/2)]

void c_solve_eigen(disort_state *ds,
                   int           lc,
                   disort_pair  *ab,
                   double       *array,
                   double       *cmu,
                   double       *cwt,
                   double       *gl,
                   int           mazim,
                   int           nn,
                   double       *ylmc,
                   double       *cc,
                   double       *evecc,
                   double       *eval,
                   double       *kk,
                   double       *gc,
                   double       *wk)
{
  int
    ier;
  register int
    iq,jq,kq,l;
  double
    alpha,beta,gpmigm,gpplgm,sum;

  /*
   * Calculate quantities in eqs. SS(5-6), STWL(8b,15,23f)
   */
  for (iq = 1; iq <= nn; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      sum = 0.;
      for (l = mazim; l <= ds->nstr-1; l++) {
        sum += GL(l,lc)*YLMC(l,iq)*YLMC(l,jq);
      }
      CC(iq,jq) = .5*sum*CWT(jq);
    }
    for (jq = 1; jq <= nn; jq++) {
      /*
       * Fill remainder of array using symmetry relations  C(-mui,muj) = C(mui,-muj) and C(-mui,-muj) = C(mui,muj)
       */
      CC(iq+nn,jq   ) = CC(iq,jq+nn);
      CC(iq+nn,jq+nn) = CC(iq,jq   );
      /*
       * Get factors of coeff. matrix of reduced eigenvalue problem
       */
      alpha      = CC(iq,jq   )/CMU(iq);
      beta       = CC(iq,jq+nn)/CMU(iq);
      AMB(iq,jq) = alpha-beta;
      APB(iq,jq) = alpha+beta;
    }
    AMB(iq,iq) -= 1./CMU(iq);
    APB(iq,iq) -= 1./CMU(iq);
  }
  /*
   * Finish calculation of coefficient matrix of reduced eigenvalue problem:
   * get matrix product (alpha+beta)*(alpha-beta); SS(12),STWL(23f)
   */
  for (iq = 1; iq <= nn; iq++) {
    for (jq = 1; jq <= nn; jq++) {
      sum = 0.;
      for (kq = 1; kq <= nn; kq++) {
        sum += APB(iq,kq)*AMB(kq,jq);
      }
      ARRAY(iq,jq) = sum;
    }
  }

  /*
   * Find (real) eigenvalues and eigenvectors
   */
  c_asymmetric_matrix(array,evecc,eval,nn,ds->nstr/2,ds->nstr,&ier,wk);

  if (ier > 0) {
    printf("\n\n asymmetric_matrix--eigenvalue no. %4d didn't converge.  Lower-numbered eigenvalues wrong.\n",ier);
    c_errmsg("asymmetric_matrix--convergence problems",DS_ERROR);
  }

  for (iq = 1; iq <= nn; iq++) {
    EVAL(iq)     = sqrt(fabs(EVAL(iq)));
    KK(iq+nn,lc) = EVAL(iq);
    /*
     * Add negative eigenvalue
     */
    KK(nn+1-iq,lc) = -EVAL(iq);
  }

  /*
   * Find eigenvectors (G+) + (G-) from SS(10) and store temporarily in APB array
   */
  for (jq = 1; jq <= nn; jq++) {
    for (iq = 1; iq <= nn; iq++) {
      sum = 0.;
      for (kq = 1; kq <= nn; kq++) {
        sum += AMB(iq,kq)*EVECC(kq,jq);
      }
      APB(iq,jq) = sum/EVAL(jq);
    }
  }
  for (jq = 1; jq <= nn; jq++) {
    for (iq = 1; iq <= nn; iq++) {
      gpplgm = APB(  iq,jq);
      gpmigm = EVECC(iq,jq);
      /*
       * Recover eigenvectors G+,G- from their sum and difference; stack them to get eigenvectors of full system
       * SS(7) (JQ = eigenvector number)
       */
      EVECC(iq,   jq) = .5*(gpplgm+gpmigm);
      EVECC(iq+nn,jq) = .5*(gpplgm-gpmigm);
      /*
       * Eigenvectors corresponding to negative eigenvalues (corresp. to reversing sign of 'k' in SS(10) )
       */
      gpplgm *= -1;
      EVECC(iq,   jq+nn)     = .5*(gpplgm+gpmigm);
      EVECC(iq+nn,jq+nn)     = .5*(gpplgm-gpmigm);
      GC(nn+iq,  nn+jq,  lc) = EVECC(iq,   jq   );
      GC(nn-iq+1,nn+jq,  lc) = EVECC(iq+nn,jq   );
      GC(nn+iq,  nn-jq+1,lc) = EVECC(iq,   jq+nn);
      GC(nn-iq+1,nn-jq+1,lc) = EVECC(iq+nn,jq+nn);
    }
  }

  return;
}

/*============================= end of c_solve_eigen() ==================*/
