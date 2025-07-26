// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_chapman_simpler() ======================*/

double c_chapman_simpler(int     lc,
                 double  taup,
                 int     nlyr,
                 double *zd,
                 double *dtau_c,
                 double  zenang,
                 double  r)
{
  register int
    id,j;
  double
    zenrad,xp,xpsinz,
    sum,fact,fact2,rj,rjp1,dhj,dsj;

  zenrad = zenang*DEG;
  xp     = r+ZD(lc)+(ZD(lc-1)-ZD(lc))*taup;
  xpsinz = xp*sin(zenrad);

  if (zenang > 90. && xpsinz < r) {
    return 1.e+20;
  }

  /*
   * Find index of layer in which the screening height lies
   */
  id = lc;
  if (zenang > 90.) {
    for (j= lc; j <= nlyr; j++) {
      if (xpsinz < (ZD(j-1)+r) && (xpsinz >= ZD(j)+r)) {
        id = j;
      }
    }
  }

  sum = 0.;
  for (j = 1; j <= id; j++) {
    fact  = 1.;
    fact2 = 1.;
    /*
     * Include factor of 2 for zenang > 90., second sum in eq. B2 (DS)
     */
    if (j > lc) {
      fact = 2.;
    }
    else if (j == lc && lc == id && zenang > 90.) {
      fact2 = -1.;
    }

    rj   = r+ZD(j-1);
    rjp1 = r+ZD(j  );
    if (j == lc && id == lc) {
      rjp1 = xp;
    }

    dhj = ZD(j-1)-ZD(j);
    if (id > lc && j == id) {
      dsj = sqrt(rj*rj-xpsinz*xpsinz);
    }
    else {
      dsj = sqrt(rj*rj-xpsinz*xpsinz)-fact2*sqrt(rjp1*rjp1-xpsinz*xpsinz);
    }
    sum += DTAU_C(j)*fact*dsj/dhj;
  }
  /*
   * Add third term in eq. B2 (DS)
   */
  if (id > lc) {
    dhj  = ZD(lc-1)-ZD(lc);
    dsj  = sqrt(xp*xp-xpsinz*xpsinz)-sqrt(SQR(ZD(lc)+r)-xpsinz*xpsinz);
    sum += DTAU_C(lc)*dsj/dhj;
  }

  return sum;
}

/*============================= end of c_chapman_simpler() ===============*/
