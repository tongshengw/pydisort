// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_asymmetric_matrix() ===================*/

/*
  Solves eigenfunction problem for real asymmetric matrix for which it
  is known a priori that the eigenvalues are real. This is an adaptation
  of a subroutine EIGRF in the IMSL library to use real instead of complex
  arithmetic, accounting for the known fact that the eigenvalues and
  eigenvectors in the discrete ordinate solution are real.

  EIGRF is based primarily on EISPACK routines.  The matrix is first
  balanced using the Parlett-Reinsch algorithm.  Then the Martin-Wilkinson
  algorithm is applied. There is a statement 'j = wk(i)' that converts a
  double precision variable to an integer variable; this seems dangerous
  to us in principle, but seems to work fine in practice.

  References:

  Dongarra, J. and C. Moler, EISPACK -- A Package for Solving Matrix
      Eigenvalue Problems, in Cowell, ed., 1984: Sources and Development of
      Mathematical Software, Prentice-Hall, Englewood Cliffs, NJ
  Parlett and Reinsch, 1969: Balancing a Matrix for Calculation of
      Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
  Wilkinson, J., 1965: The Algebraic Eigenvalue Problem, Clarendon Press,
      Oxford

   I N P U T    V A R I A B L E S:

       aa    :  input asymmetric matrix, destroyed after solved
        m    :  order of aa
       ia    :  first dimension of aa
    ievec    :  first dimension of evec

   O U T P U T    V A R I A B L E S:

       evec  :  (unnormalized) eigenvectors of aa (column j corresponds to EVAL(J))
       eval  :  (unordered) eigenvalues of aa (dimension m)
       ier   :  if != 0, signals that EVAL(ier) failed to converge;
                   in that case eigenvalues ier+1,ier+2,...,m  are
                   correct but eigenvalues 1,...,ier are set to zero.

   S C R A T C H   V A R I A B L E S:

       wk    :  work area (dimension at least 2*m)

   Called by- c_solve_eigen
   Calls- c_errmsg
 -------------------------------------------------------------------*/

void c_asymmetric_matrix(double *aa,
                         double *evec,
                         double *eval,
                         int     m,
                         int     ia,
                         int     ievec,
                         int    *ier,
                         double *wk)
{
  const double
   c1 =    .4375,
   c2 =    .5,
   c3 =    .75,
   c4 =    .95,
   c5 =  16.,
   c6 = 256.;
  int
    noconv,notlas,
    i,ii,in,j,k,ka,kkk,l,lb=0,lll,n,n1,n2;
  double
    col,discri,f,g,h,p=0,q=0,r=0,repl,rnorm,row,
    s,scale,sgn,t,tol,uu,vv,w,x,y,z;

  *ier = 0;
  tol = DBL_EPSILON;
  if (m < 1 || ia < m || ievec < m) {
    c_errmsg("asymmetric_matrix--bad input variable(s)",DS_ERROR);
  }

  /*
   * Handle 1x1 and 2x2 special cases
   */
  if (m == 1) {
    EVAL(1)   = AA(1,1);
    EVEC(1,1) = 1.;
    return;
  }
  else if (m == 2) {
    discri = SQR(AA(1,1)-AA(2,2))+4.*AA(1,2)*AA(2,1);
    if(discri < 0.) {
      c_errmsg("asymmetric_matrix--complex evals in 2x2 case",DS_ERROR);
    }
    sgn = 1.;
    if (AA(1,1) < AA(2,2)) {
     sgn = -1.;
    }
    EVAL(1)   = .5*(AA(1,1)+AA(2,2)+sgn*sqrt(discri));
    EVAL(2)   = .5*(AA(1,1)+AA(2,2)-sgn*sqrt(discri));
    EVEC(1,1) = 1.;
    EVEC(2,2) = 1.;
    if (AA(1,1) == AA(2,2) && (AA(2,1) == 0. || AA(1,2) == 0.)) {
      rnorm     = fabs(AA(1,1))+fabs(AA(1,2))+fabs(AA(2,1))+fabs(AA(2,2));
      w         = tol*rnorm;
      EVEC(2,1) =  AA(2,1)/w;
      EVEC(1,2) = -AA(1,2)/w;
    }
    else {
      EVEC(2,1) = AA(2,1)/(EVAL(1)-AA(2,2));
      EVEC(1,2) = AA(1,2)/(EVAL(2)-AA(1,1));
    }
    return;
  }

  /*
   * Initialize output variables
   */
  *ier = 0;
  memset(eval,0,m*sizeof(double));
  memset(evec,0,ievec*ievec*sizeof(double));
  for (i = 1; i <= m; i++) {
    EVEC(i,i) = 1.;
  }

  /*
   * Balance the input matrix and reduce its norm by diagonal similarity transformation stored in wk;
   * then search for rows isolating an eigenvalue and push them down.
   */
  rnorm = 0.;
  l     = 1;
  k     = m;

S50:

  kkk = k;
  for (j = kkk; j >= 1; j--) {
    row = 0.;
    for (i = 1; i <= k; i++) {
      if (i != j) {
        row += fabs(AA(j,i));
      }
    }
    if (row == 0.) {
      WK(k) = (double)j;
      if (j != k) {
        for (i = 1; i <= k; i++) {
          repl    = AA(i,j);
          AA(i,j) = AA(i,k);
          AA(i,k) = repl;
        }
        for (i = l; i <= m; i++) {
          repl    = AA(j,i);
          AA(j,i) = AA(k,i);
          AA(k,i) = repl;
        }
      }
      k--;
      goto S50;
    }
  }

  /*
   * Search for columns isolating an eigenvalue and push them left.
   */

S100:

  lll = l;
  for (j = lll; j <= k; j++) {
    col = 0.;
    for (i = l; i <= k; i++) {
      if (i != j) {
        col += fabs(AA(i,j));
      }
    }
    if (col == 0.) {
      WK(l) = (double)j;
      if (j != l) {
        for (i = 1; i <= k; i++) {
          repl    = AA(i,j);
          AA(i,j) = AA(i,l);
          AA(i,l) = repl;
        }
        for (i = l; i <= m; i++) {
          repl    = AA(j,i);
          AA(j,i) = AA(l,i);
          AA(l,i) = repl;
        }
      }
      l++;
      goto S100;
    }
  }

  /*
   * Balance the submatrix in rows L through K
   */
  for (i = l; i <= k; i++) {
    WK(i) = 1.;
  }

  noconv = TRUE;
  while (noconv) {
    noconv = FALSE;
    for (i = l; i <= k; i++) {
      col = 0.;
      row = 0.;
      for (j = l; j <= k; j++) {
        if (j != i) {
          col += fabs(AA(j,i));
          row += fabs(AA(i,j));
        }
      }

      f = 1.;
      g = row/c5;
      h = col+row;

      while (col < g) {
        f   *= c5;
        col *= c6;
      }

      g = row*c5;

      while (col >= g) {
        f   /= c5;
        col /= c6;
      }

      /*
       * Now balance
       */
      if ((col+row)/f < c4*h) {
        WK(i)  *= f;
        noconv  = TRUE;
        for (j = l; j <= m; j++) {
          AA(i,j) /= f;
        }
        for (j = 1; j <= k; j++) {
          AA(j,i) *= f;
        }
      }
    }
  }

  if (k-1 >= l+1) {
    /*
     * Transfer A to a Hessenberg form.
     */
    for (n = l+1; n <= k-1; n++) {
      h       = 0.;
      WK(n+m) = 0.;
      scale   = 0.;
      /*
       * Scale column
       */
      for (i = n; i <= k; i++) {
        scale += fabs(AA(i,n-1));
      }
      if (scale != 0.) {
        for (i = k; i >= n; i--) {
          WK(i+m)  = AA(i,n-1)/scale;
          h       += SQR(WK(i+m));
        }
        g        = -F77_SIGN(sqrt(h),WK(n+m));
        h       -= WK(n+m)*g;
        WK(n+m) -= g;
        /*
         * Form (I-(U*UT)/H)*A
         */
        for (j = n; j <= m; j++) {
          f = 0.;
          for (i = k; i >= n; i--) {
            f += WK(i+m)*AA(i,j);
          }
          for (i = n; i <= k; i++) {
            AA(i,j) -= WK(i+m)*f/h;
          }
        }
        /*
         * Form (i-(u*ut)/h)*a*(i-(u*ut)/h)
         */
        for (i = 1; i <= k; i++) {
          f = 0.;
          for (j = k; j >= n; j--) {
            f += WK(j+m)*AA(i,j);
          }
          for (j = n; j <= k; j++) {
            AA(i,j) -= WK(j+m)*f/h;
          }
        }
        WK(n+m)   *= scale;
        AA(n,n-1)  = scale*g;
      }
    }

    for (n = k-2; n >= l; n--) {
      n1 = n+1;
      n2 = n+2;
      f = AA(n+1,n);
      if( f != 0.) {
        f *= WK(n+1+m);
        for (i = n+2; i <= k; i++) {
          WK(i+m) = AA(i,n);
        }
        if (n+1 <= k) {
          for (j = 1; j <= m; j++) {
            g = 0.;
            for (i = n+1; i <= k; i++) {
              g += WK(i+m)*EVEC(i,j);
            }
            g /= f;
            for (i = n+1; i <= k; i++) {
              EVEC(i,j) += g*WK(i+m);
            }
          }
        }
      }
    }
  }

  n = 1;
  for (i = 1; i <= m; i++) {
    for (j = n; j <= m; j++) {
      rnorm += fabs(AA(i,j));
    }
    n = i;
    if (i < l || i > k) {
      EVAL(i) = AA(i,i);
    }
  }

  n = k;
  t = 0.;
  /*
   * Search for next eigenvalues
   */

S400:

  if (n < l) {
    goto S550;
  }

  in = 0;
  n1 = n-1;
  n2 = n-2;

  /*
   * Look for single small sub-diagonal element
   */

S410:

  for (i = l; i <= n; i++) {
    lb = n+l-i;
    if (lb == l) {
      break;
    }
    s = fabs(AA(lb-1,lb-1))+fabs(AA(lb,lb));
    if (s == 0.) {
      s = rnorm;
    }
    if (fabs(AA(lb,lb-1)) <= tol*s) {
      break;
    }
  }

  x = AA(n,n);
  if (lb == n) {
    /*
     * One eigenvalue found
     */
    AA(n,n) = x+t;
    EVAL(n) = AA(n,n);
    n       = n1;
    goto S400;
  }

  y = AA(n1,n1);
  w = AA(n,n1)*AA(n1,n);

  if (lb == n1) {
    /*
     * Two eigenvalues found
     */
    p         = (y-x)*c2;
    q         = p*p+w;
    z         = sqrt(fabs(q));
    AA(n,n)   = x+t;
    x         = AA(n,n);
    AA(n1,n1) = y+t;
    /*
     * Real pair
     */
    z        = p+F77_SIGN(z,p);
    EVAL(n1) = x+z;
    EVAL(n)  = EVAL(n1);

    if (z != 0.) {
      EVAL(n) = x-w/z;
    }
    x = AA(n,n1);
    /*
     * Employ scale factor in case X and Z are very small
     */
    r = sqrt(x*x+z*z);
    p = x/r;
    q = z/r;
    /*
     * Row modification
     */
    for (j = n1; j <= m; j++) {
      z        = AA(n1,j);
      AA(n1,j) =  q*z+p*AA(n,j);
      AA(n, j) = -p*z+q*AA(n,j);
    }
    /*
     * Column modification
     */
    for (i = 1; i <= n; i++) {
      z        = AA(i,n1);
      AA(i,n1) =  q*z+p*AA(i,n);
      AA(i,n ) = -p*z+q*AA(i,n);
    }
    /*
     * Accumulate transformations
     */
    for (i = l; i <= k; i++) {
      z          = EVEC(i,n1);
      EVEC(i,n1) =  q*z+p*EVEC(i,n);
      EVEC(i,n ) = -p*z+q*EVEC(i,n);
    }
    n = n2;
    goto S400;
  }

  if (in == 30) {
    /*
     * No convergence after 30 iterations; set error indicator to
     * the index of the current eigenvalue, and return.
     */
    *ier = n;
    return;
  }

  /*
   * Form shift
   */
  if (in == 10 || in == 20) {
    t += x;
    for (i = l; i <= n; i++) {
      AA(i,i) -= x;
    }
    s = fabs(AA(n,n1))+fabs(AA(n1,n2));
    x = c3*s;
    y = x;
    w = -c1*s*s;
  }

  in++;

  /*
   * Look for two consecutive small sub-diagonal elements
   */
  for (j = lb; j <= n2; j++) {
    i  = n2+lb-j;
    z  = AA(i,i);
    r  = x-z;
    s  = y-z;
    p  = (r*s-w)/AA(i+1,i)+AA(i,i+1);
    q  = AA(i+1,i+1)-z-r-s;
    r  = AA(i+2,i+1);
    s  = fabs(p)+fabs(q)+fabs(r);
    p /= s;
    q /= s;
    r /= s;

    if (i == lb) {
      break;
    }

    uu = fabs(AA(i,i-1))*(fabs(q)+fabs(r));
    vv = fabs(p)*(fabs(AA(i-1,i-1))+fabs(z)+fabs(AA(i+1,i+1)));

    if (uu <= tol*vv) {
      break;
    }
  }

  AA(i+2,i) = 0.;
  for (j = i+3; j <= n; j++) {
    AA(j,j-2) = 0.;
    AA(j,j-3) = 0.;
  }

  /*
   * Double QR step involving rows K to N and columns M to N
   */
  for (ka = i; ka <= n1; ka++) {
    notlas = (ka != n1);
    if (ka == i) {
      s = F77_SIGN(sqrt(p*p+q*q+r*r),p);
      if (lb != i) {
        AA(ka,ka-1) *= -1;
      }
    }
    else {
      p = AA(ka,  ka-1);
      q = AA(ka+1,ka-1);
      r = 0.;
      if (notlas) {
        r = AA(ka+2,ka-1);
      }
      x = fabs(p)+fabs(q)+fabs(r);
      if (x == 0.) {
        continue;
      }
      p /= x;
      q /= x;
      r /= x;
      s  = F77_SIGN(sqrt(p*p+q*q+r*r),p);

      AA(ka,ka-1) = -s*x;
    }

    p += s;
    x  = p/s;
    y  = q/s;
    z  = r/s;
    q /= p;
    r /= p;

    /*
     * Row modification
     */
    for (j = ka; j <= m; j++) {
      p = AA(ka,j)+q*AA(ka+1,j);
      if (notlas) {
        p          += r*AA(ka+2,j);
        AA(ka+2,j) -= p*z;
      }
      AA(ka+1,j) -= p*y;
      AA(ka,  j) -= p*x;
    }

    /*
     * Column modification
     */
    for (ii = 1; ii <= IMIN(n,ka+3); ii++) {
      p = x*AA(ii,ka)+y*AA(ii,ka+1);
      if (notlas) {
        p           += z*AA(ii,ka+2);
        AA(ii,ka+2) -= p*r;
      }
      AA(ii,ka+1) -= p*q;
      AA(ii,ka  ) -= p;
    }

    /*
     * Accumulate transformations
     */
    for (ii = l; ii <= k; ii++) {
      p = x*EVEC(ii,ka)+y*EVEC(ii,ka+1);
      if (notlas) {
        p             += z*EVEC(ii,ka+2);
        EVEC(ii,ka+2) -= p*r;
      }
      EVEC(ii,ka+1) -= p*q;
      EVEC(ii,ka  ) -= p;
    }
  }

  goto S410;

  /*
   * All evals found, now backsubstitute real vector
   */

S550:

  if (rnorm != 0.) {
    for (n = m; n >= 1; n--) {
      n2      = n;
      AA(n,n) = 1.;
      for (i = n-1; i >= 1; i--) {
        w = AA(i,i)-EVAL(n);
        if (w == 0.) {
          w = tol*rnorm;
        }
        r = AA(i,n);
        for (j = n2; j <= n-1; j++) {
          r += AA(i,j)*AA(j,n);
        }
        AA(i,n) = -r/w;
        n2      = i;
      }
    }
    /*
     * End backsubstitution vectors of isolated evals
     */
    for (i = 1; i <= m; i++) {
      if (i < l || i > k) {
        for (j = i; j <= m; j++) {
          EVEC(i,j) = AA(i,j);
        }
      }
    }
    /*
     * Multiply by transformation matrix
     */
    if (k != 0) {
      for (j = m; j >= l; j--) {
        for (i = l; i <= k; i++) {
          z = 0.;
          for (n = l; n <= IMIN(j,k); n++) {
            z += EVEC(i,n)*AA(n,j);
          }
          EVEC(i,j) = z;
        }
      }
    }
  }
  for (i = l; i <= k; i++) {
    for (j = 1; j <= m; j++) {
      EVEC(i,j) *= WK(i);
    }
  }

  /*
   * Interchange rows if permutations occurred
   */
  for (i = l-1; i >= 1; i--) {
    j = WK(i);
    if (i != j) {
      for (n = 1; n <= m; n++) {
        repl      = EVEC(i,n);
        EVEC(i,n) = EVEC(j,n);
        EVEC(j,n) = repl;
      }
    }
  }
  for (i = k+1; i <= m; i++) {
    j = WK(i);
    if (i != j) {
      for (n = 1; n <= m; n++) {
        repl      = EVEC(i,n);
        EVEC(i,n) = EVEC(j,n);
        EVEC(j,n) = repl;
      }
    }
  }

  return;
}

/*============================= end of c_asymmetric_matrix() ============*/
