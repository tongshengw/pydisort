// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_chapman() ==============================*/

/*
 Calculates the Chapman factor.

 I n p u t       v a r i a b l e s:

      lc        : Computational layer
      taup      :
      tauc      :
      nlyr      : Number of layers in atmospheric model
      zd(lc)    : lc = 0, nlyr. zd(lc) is distance from bottom
                  surface to top of layer lc. zd(nlyr) = 0. km
      dtau_c    : Optical thickness of layer lc (un-delta-m-scaled)
      zenang    : Solar zenith angle as seen from bottom surface
      r         : Radial parameter, see Velinow & Kostov (2001). NOTE: Use the same dimension as zd,
                  for instance both in km.

 O u t p u t      v a r i a b l e s:

      ch        : Chapman-factor. In a pseudo-spherical atmosphere, replace exp(-tau/umu0) by exp(-ch(lc)) in the
                  beam source in

 I n t e r n a l     v a r i a b l e s:

      dhj       : delta-h-sub-j in eq. B2 (DS)
      dsj       : delta-s-sub-j in eq. B2 (DS)
      fact      : =1 for first  sum in eq. B2 (DS)
                  =2 for second sum in eq. B2 (DS)
      rj        : r-sub-j   in eq. B1 (DS)
      rjp1      : r-sub-j+1 in eq. B1 (DS)
      xpsinz    : The length of the line OG in Fig. 1, (DS)


 NOTE: Assumes a spherical planet. One might consider generalizing following
       Velinow YPI, Kostov VI, 2001, Generalization on Chapman Function for the Atmosphere of an Oblate Rotating Planet,
         Comptes Rendus de l'Academie Bulgare des Sciences 54, 29-34.
*/

double c_chapman(int     lc,
                 double  taup,
                 double *tauc,
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

/*============================= end of c_chapman() =======================*/
