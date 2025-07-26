// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_solns() =========================*/

/*
    Calculates the homogenous and particular solutions to the
    radiative transfer equation in the two-stream approximation,
    for each layer in the medium.

    I n p u t     v a r i a b l e s:

      ds         : 'Disort' state variables
      ch         : Chapman correction factor
      chtau      :
      cmu        : Abscissa for gauss quadrature over angle cosine
      ncut       : Number of computational layer where absorption optical depth exceeds -abscut-
      oprim      : Delta-m scaled single scattering albedo
      pkag,c     : Planck function in each layer
      flag.spher : spher = true => spherical geometry invoked
      taucpr     : Cumulative optical depth (delta-m-scaled)
      ggprim     :

   O u t p u t     v a r i a b l e s:

      kk         :  Eigenvalues
      rr         :  Eigenvectors at polar quadrature angles
      ts         :  twostr_xyz structure variables (see cdisort.h)
  ----------------------------------------------------------------------*/

void c_twostr_solns(disort_state *ds,
                    double       *ch,
                    double       *chtau,
                    double        cmu,
                    int           ncut,
                    double       *oprim,
                    double       *pkag,
                    double       *pkagc,
                    double       *taucpr,
                    double       *ggprim,
                    double       *kk,
                    double       *rr,
                    twostr_xyz   *ts)
{
  register int
    lc;
  static int
    initialized = FALSE;
  static double
    big,large,small,little;
  double
    q_1,q_2,qq,q0a,q0,q1a,q2a,q1,q2,
    deltat,denomb,z0p,z0m,arg,sgn,fact3,denomp,
    beta,fact1,fact2;

  if (!initialized) {
    /*
     * The calculation of the particular solutions require some care; small,little,
       big, and large have been set so that no problems should occur in double precision.
     */
    small  = 1.e+30*DBL_MIN;
    little = 1.e+20*DBL_MIN;
    big    = sqrt(DBL_MAX)/1.e+10;
    large  = log(DBL_MAX)-20.;

    initialized = TRUE;
  }

  /*----------------  Begin loop on computational layers  ---------------------*/

  for (lc = 1; lc <= ncut; lc++) {
    /*
     * Calculate eigenvalues -kk- and eigenvector -rr-, eqs. KST(20-21)
     */
    beta   = 0.5*(1.-3.*GGPRIM(lc)*cmu*cmu);
    fact1  = 1.-OPRIM(lc);
    fact2  = 1.-OPRIM(lc)+2.*OPRIM(lc)*beta;
    KK(lc) = (1./cmu)*sqrt(fact1*fact2);
    RR(lc) = (sqrt(fact2)-sqrt(fact1))/(sqrt(fact2)+sqrt(fact1));

    if (ds->bc.fbeam > 0.) {
      /*
       * Set coefficients in KST(22) for beam source
       */
      q_1 = ds->bc.fbeam/(4.*M_PI)*OPRIM(lc)*(1.-3.*GGPRIM(lc)*cmu*ds->bc.umu0);
      q_2 = ds->bc.fbeam/(4.*M_PI)*OPRIM(lc)*(1.+3.*GGPRIM(lc)*cmu*ds->bc.umu0);

      if (ds->bc.umu0 >= 0.) {
        qq = q_2;
      }
      else {
        qq = q_1;
      }

      if (ds->flag.spher) {
        q0a = exp(-CHTAU(lc-1));
        q0  = q0a*qq;
        if (q0 <= small) {
          q1a = 0.;
          q2a = 0.;
        }
        else {
          q1a = exp(-CHTAU(lc-1  ));
          q2a = exp(-CHTAU(lc));
        }
      }
      else {
        q0a = exp(-TAUCPR(lc-1)/ds->bc.umu0);
        q0  = q0a*qq;
        if (q0 <= small) {
          q1a = 0.;
          q2a = 0.;
        }
        else {
          q1a = exp(-(TAUCPR(lc-1)+TAUCPR(lc))/(2.*ds->bc.umu0));
          q2a = exp(-TAUCPR(lc)/ds->bc.umu0);
        }
      }
      q1 = q1a*qq;
      q2 = q2a*qq;

      /*
       * Calculate alpha coefficient
       */
      deltat     = TAUCPR(lc)-TAUCPR(lc-1);
      ZB_A(lc)   = 1./CH(lc);
      if (fabs(ZB_A(lc)*TAUCPR(lc-1)) > large || fabs(ZB_A(lc)*TAUCPR(lc)) > large) {
        ZB_A(lc) = 0.;
      }

      /*
       * Dither alpha if it is close to an eigenvalue
       */
      denomb = fact1*fact2-SQR(ZB_A(lc)*cmu);
      if (denomb < 1.e-03) {
        ZB_A(lc) = 1.02*ZB_A(lc);
      }
      q0 = q0a*q_1;
      q2 = q2a*q_1;

      /*
       * Set constants in eq. KST(22)
       */
      if (deltat < 1.e-07) {
        XB_1D(lc) = 0.;
      }
      else {
        XB_1D(lc) = 1./deltat*(q2*exp(ZB_A(lc)*TAUCPR(lc))-q0*exp(ZB_A(lc)*TAUCPR(lc-1)));
      }
      XB_0D(lc) = q0*exp(ZB_A(lc)*TAUCPR(lc-1))-XB_1D(lc)*TAUCPR(lc-1);
      q0        = q0a*q_2;
      q2        = q2a*q_2;

      if (deltat < 1.e-07) {
        XB_1U(lc) = 0.;
      }
      else {
        XB_1U(lc) = 1./deltat*(q2*exp(ZB_A(lc)*TAUCPR(lc))-q0*exp(ZB_A(lc)*TAUCPR(lc-1)));
      }
      XB_0U(lc) = q0*exp(ZB_A(lc)*TAUCPR(lc-1))-XB_1U(lc)*TAUCPR(lc-1);

      /*
       * Calculate particular solutions for incident beam source in pseudo-spherical geometry, eqs. KST(24-25)
       */
      denomb    = fact1*fact2-SQR(ZB_A(lc)*cmu);
      YB_1D(lc) = (OPRIM(lc)*beta*XB_1D(lc)+(1.-OPRIM(lc)*(1.-beta)+ZB_A(lc)*cmu)*XB_1U(lc))/denomb;
      YB_1U(lc) = (OPRIM(lc)*beta*XB_1U(lc)+(1.-OPRIM(lc)*(1.-beta)-ZB_A(lc)*cmu)*XB_1D(lc))/denomb;
      z0p       = XB_0U(lc)-cmu*YB_1D(lc);
      z0m       = XB_0D(lc)+cmu*YB_1U(lc);
      YB_0D(lc) = (OPRIM(lc)*beta*z0m+(1.-OPRIM(lc)*(1.-beta)+ZB_A(lc)*cmu)*z0p)/denomb;
      YB_0U(lc) = (OPRIM(lc)*beta*z0p+(1.-OPRIM(lc)*(1.-beta)-ZB_A(lc)*cmu)*z0m)/denomb;
    }

    if(ds->flag.planck) {
      /*
       * Set coefficients in KST(22) for thermal source
       * Calculate alpha coefficient
       */
      q0     = (1.-OPRIM(lc))*PKAG(lc-1);
      q1     = (1.-OPRIM(lc))*PKAGC(lc);
      q2     = (1.-OPRIM(lc))*PKAG(lc);
      deltat = TAUCPR(lc)-TAUCPR(lc-1);

      if ((q2 < q0*1.e-02 || q2 <= little) && q1 > little && q0 > little) {
        /*
         * Case 1: source small at bottom layer; alpha eq. KS(50)
         */
        ZP_A(lc) = MIN(2./deltat*log(q0/q1),big);
        if (ZP_A(lc)*TAUCPR(lc-1) >= log(big)) {
          XP_0(lc) = big;
        }
        else {
          XP_0(lc) = q0;
        }
        XP_1(lc) = 0.;
      }
      else if ((q2 <= q1*1.e-02 || q2 <= little) && (q1 <= q0*1.e-02 || q1 <= little) && q0 > little) {
        /*
         * Case 2: Source small at center and bottom of layer
         */
        ZP_A(lc) = big/TAUCPR(ncut);
        XP_0(lc) = q0;
        XP_1(lc) = 0.;
      }
      else if (q2 <= little && q1 <= little && q0 <= little) {
        /*
         * Case 3: All sources zero
         */
        ZP_A(lc) = 0.;
        XP_0(lc) = 0.;
        XP_1(lc) = 0.;
      }
      else if ( ( fabs((q2-q0)/q2) < 1.e-04 && fabs((q2-q1)/q2) < 1.e-04 ) || deltat < 1.e-04) {
        /*
         * Case 4: Sources same at center, bottom and top of layer or layer optically very thin
         */
        ZP_A(lc) = 0.;
        XP_0(lc) = q0;
        XP_1(lc) = 0.;
      }
      else {
        /*
         *  Case 5: Normal case
         */
        arg = MAX(SQR(q1/q2)-q0/q2,0.);
        /*
         * alpha eq. (44). For source that has its maximum at the top of the layer, use negative solution
         */
        sgn = 1.;
        if (PKAG(lc-1) > PKAG(lc)) {
         sgn = -1.;
        }
        fact3 = log(q1/q2+sgn*sqrt(arg));

        /* Be careful with log of numbers close to one */
        if (fabs(fact3) <= 0.005) {
          /* numbers close to one */
          q1    = 0.99*q1;
          fact3 = log(q1/q2+sgn*sqrt(arg));
        }

        ZP_A(lc) = 2./deltat*fact3;
        if (fabs(ZP_A(lc)*TAUCPR(lc)) > log(DBL_MAX)-log(q0*100.)) {
          ZP_A(lc) = 0.;
        }

        /*
         * Dither alpha if it is close to an eigenvalue
         */
        denomp = fact1*fact2-SQR(ZP_A(lc)*cmu);
        if (denomp < 1.e-03) {
          ZP_A(lc) *= 1.01;
        }

        /*
         * Set constants in eqs. KST(22)
         */
        if(deltat < 1.e-07) {
          XP_1(lc) = 0.;
        }
        else {
          XP_1(lc) = 1./deltat*(q2*exp(ZP_A(lc)*TAUCPR(lc))-q0*exp(ZP_A(lc)*TAUCPR(lc-1)));
        }
        XP_0(lc) = q0*exp(ZP_A(lc)*TAUCPR(lc-1))-XP_1(lc)*TAUCPR(lc-1);
      }

      /*
       * Calculate particular solutions eqs. KST(24-25) for internal thermal so
       */
      denomp    = fact1*fact2-SQR(ZP_A(lc)*cmu);
      YP_1D(lc) = (OPRIM(lc)*beta*XP_1(lc)+(1.-OPRIM(lc)*(1.-beta)+ZP_A(lc)*cmu)*XP_1(lc))/denomp;
      YP_1U(lc) = (OPRIM(lc)*beta*XP_1(lc)+(1.-OPRIM(lc)*(1.-beta)-ZP_A(lc)*cmu)*XP_1(lc))/denomp;
      z0p       = XP_0(lc)-cmu*YP_1D(lc);
      z0m       = XP_0(lc)+cmu*YP_1U(lc);
      YP_0D(lc) = (OPRIM(lc)*beta*z0m+(1.-OPRIM(lc)*(1.-beta)+ZP_A(lc)*cmu)*z0p)/denomp;
      YP_0U(lc) = (OPRIM(lc)*beta*z0p+(1.-OPRIM(lc)*(1.-beta)-ZP_A(lc)*cmu)*z0m)/denomp;
    }
  }

  return;
}

/*============================= end of c_twostr_solns() ==================*/
