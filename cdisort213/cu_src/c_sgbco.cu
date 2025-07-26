// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_sgbco() ================================*/

/*
     Factors a real band matrix by Gaussian elimination and estimates the
     condition of the matrix.
     Revision date:  8/1/82
     Author:  Moler, C.B. (Univ. of New Mexico)

     If  RCOND  is not needed, sgbfa is slightly faster.
     To solve  A*X = B , follow sgbco by sgbsl.

     Inputs:
        abd     double(LDA,N), contains the matrix in band storage.
                The columns of the matrix are stored in the columns of abd
                and the diagonals of the matrix are stored in rows
                ml+1 through 2*ml+mu+1 of  abd.
                See the comments below for details.
        lda     int, the leading dimension of the array abd.
                lda must be >= 2*ml+mu+1.
        n       int,the order of the original matrix.
        ml      int, number of diagonals below the main diagonal.
                0 <= ml < n.
        mu      int, number of diagonals above the main diagonal.
                0 <= mu < n.
                more efficient if  ml <= mu.

     Outputs:
        abd     an upper triangular matrix in band storage and
                the multipliers which were used to obtain it.
                The factorization can be written  A = L*U  where
                L  is a product of permutation and unit lower
                triangular matrices and  U  is upper triangular.
        ipvt    int[n], an integer vector of pivot indices.
        rcond   double, an estimate of the reciprocal condition of A.
                For the system  A*X = B, relative perturbations
                in A and B of size epsilon may cause relative
                perturbations in  X  of size  epsilon/rcond.
                If rcond  is so small that the logical expression
                   1.+RCOND == 1.
                is true, then  A  may be singular to working
                precision.  In particular, rcond is zero if exact
                singularity is detected or the estimate underflows.
        z       double[n], a work vector whose contents are usually
                unimportant. If A is close to a singular matrix, then
                z is an approximate null vector in the sense that
                norm(a*z) = rcond*norm(a)*norm(z).

     Band storage:
           If A is a band matrix, the following program segment
           will set up the input (with unit-offset arrays):
                   ml = (band width below the diagonal)
                   mu = (band width above the diagonal)
                   m = ml+mu+1
                   for (j = 1; j <= n; j++) {
                     i1 = IMAX(1,j-mu);
                     i2 = IMIN(n,j+ml);
                     for (i = i1; i <= i2; i++) {
                       k = i-j+m;
                       ABD(K,J) = A(I,J);
                     }
                   }
           This uses rows ml+1 through 2*ml+mu+1 of abd.
           In addition, the first ml rows in abd are used for
           elements generated during the triangularization.
           The total number of rows needed in abd is 2*ml+mu+1.
           The ml+mu by ml+mu upper left triangle and the
           ml by ml lower right triangle are not referenced.

     Example:  if the original matrix is

           11 12 13  0  0  0
           21 22 23 24  0  0
            0 32 33 34 35  0
            0  0 43 44 45 46
            0  0  0 54 55 56
            0  0  0  0 65 66

      then  n = 6, ml = 1, mu = 2, lda >= 5  and abd should contain
            *  *  *  +  +  +  , * = not used
            *  * 13 24 35 46  , + = used for pivoting
            * 12 23 34 45 56
           11 22 33 44 55 66
           21 32 43 54 65  *

 --------------------------------------------------------------------*/

void c_sgbco(double *abd,
             int     lda,
             int     n,
             int     ml,
             int     mu,
             int    *ipvt,
             double *rcond,
             double *z)
{
  int
    info;
  register int
    is,j,ju,k,kb,kp1,l,la,lm,lz,m,mm;
  double
    anorm,ek,s,sm,t,wk,wkm,ynorm;

  /*
   * compute 1-norm of A
   */
  anorm = 0.;
  l  = ml+1;
  is = l+mu;
  for (j = 1; j <= n; j++) {
    anorm = MAX(anorm,c_sasum(l,&ABD(is,j)));
    if (is > ml+1) {
      is--;
    }
    if (j <= mu) {
      l++;
    }
    if (j >= n-ml) {
      l--;
    }
  }
  /*
   * factor
   */
  c_sgbfa(abd,lda,n,ml,mu,ipvt,&info);
  /*
   * rcond = 1/(norm(A)*(estimate of norm(inverse(A)))) .
   * estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
   * trans(A) is the transpose of A.  The components of E are
   * chosen to cause maximum local growth in the elements of W where
   * trans(U)*W = E. The vectors are frequently rescaled to avoid overflow.
   * solve trans(U)*W = E
   */
  ek = 1.;

  memset(z,0,n*sizeof(double));

  m  = ml+mu+1;
  ju = 0;
  for (k = 1; k <= n; k++) {
    if (Z(k) != 0.) {
      ek = F77_SIGN(ek,-Z(k));
    }
    if (fabs(ek-Z(k)) > fabs(ABD(m,k))) {
      s = fabs(ABD(m,k))/fabs(ek-Z(k));
      c_sscal(n,s,z);
      ek *= s;
    }
    wk  =  ek-Z(k);
    wkm = -ek-Z(k);
    s   = fabs(wk);
    sm  = fabs(wkm);
    if (ABD(m,k) != 0.) {
      wk  /= ABD(m,k);
      wkm /= ABD(m,k);
    }
    else {
      wk  = 1.;
      wkm = 1.;
    }
    kp1 = k+1;
    ju  = IMIN(IMAX(ju,mu+IPVT(k)),n);
    mm  = m;
    if (kp1 <= ju) {
      for (j = kp1; j <= ju; j++) {
        mm--;
        sm   += fabs(Z(j)+wkm*ABD(mm,j));
        Z(j) += wk*ABD(mm,j);
        s    += fabs(Z(j));
      }
      if (s < sm) {
        t  = wkm-wk;
        wk = wkm;
        mm = m;
        for (j = kp1; j <= ju; j++) {
          mm--;
          Z(j) += t*ABD(mm,j);
        }
      }
    }
    Z(k) = wk;
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  /*
   * solve trans(L)*Y = W
   */
  for (kb = 1; kb <= n; kb++) {
    k  = n+1-kb;
    lm = IMIN(ml,n-k);
    if (k < n) {
      Z(k) += c_sdot(lm,&ABD(m+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
    }

    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  ynorm = 1.;
  /*
   * solve L*V = Y
   */
  for (k = 1; k <= n; k++) {
    l    = IPVT(k);
    t    = Z(l);
    Z(l) = Z(k);
    Z(k) = t;
    lm   = IMIN(ml,n-k);
    if (k < n) {
      c_saxpy(lm,t,&ABD(m+1,k),&Z(k+1));
    }
    if (fabs(Z(k)) > 1.) {
      s = 1./fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
  }

  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  ynorm *= s;
  /*
   * solve  U*Z = W
   */
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (fabs(Z(k)) > fabs(ABD(m,k))) {
      s = fabs(ABD(m,k))/fabs(Z(k));
      c_sscal(n,s,z);
      ynorm *= s;
    }
    if (ABD(m,k) != 0.) {
      Z(k) /= ABD(m,k);
    }
    else {
      Z(k) = 1.;
    }
    lm = IMIN(k,m)-1;
    la = m-lm;
    lz = k-lm;
    t  = -z[k-1];
    c_saxpy(lm,t,&ABD(la,k),&Z(lz));
  }

  /*
   * make znorm = 1.
   */
  s = 1./c_sasum(n,z);
  c_sscal(n,s,z);

  ynorm *= s;
  if(anorm != 0.) {
    *rcond = ynorm/anorm;
  }
  else {
    *rcond = 0.;
  }

  return;
}

/*============================= end of c_sgbco() =========================*/
