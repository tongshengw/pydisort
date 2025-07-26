// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_user_intensities() ====================*/

/*
   Computes intensity components at user output angles for azimuthal
   expansion terms in eq. SD(2), STWL(6)

   I N P U T    V A R I A B L E S:

       ds     :  Disort state variables
       bplanck:  Integrated Planck function for emission from
                 bottom boundary
       cmu    :  Abscissae for Gauss quadrature over angle cosine
       cwt    :  Weights for Gauss quadrature over angle cosine
       delm0  :  Kronecker delta, delta-sub-M0
       emu    :  Surface directional emissivity (user angles)
       expbea :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
       gc     :  Eigenvectors at polar quadrature angles, SC(1)
       gu     :  Eigenvectors interpolated to user polar angles
                    (i.e., G in eq. SC(1) )
       kk     :  Eigenvalues of coeff. matrix in eq. SS(7), STWL(23b)
       layru  :  Layer number of user level UTAU
       ll     :  Constants of integration in eq. SC(1), obtained
                 by solving scaled version of eq. SC(5);
                 exponential term of eq. SC(12) not included
       lyrcut :  Logical flag for truncation of computational layer
       mazim  :  Order of azimuthal component
       ncut   :  Total number of computational layers considered
       nn     :  Order of double-Gauss quadrature (NSTR/2)
       rmu    :  Surface bidirectional reflectivity (user angles)
       taucpr :  Cumulative optical depth (delta-M-Scaled)
       tplanck:  Integrated Planck function for emission from
                 top boundary
       utaupr :  Optical depths of user output levels in delta-M
                 coordinates;  equal to UTAU if no delta-M
       zgu    :  General source function at user angles
       zu     :  Z-sub-zero, Z-sub-one in eq. SS(16) interpolated to user angles from an equation derived from SS(16),
                 Y-sub-zero, Y-sub-one on STWL(26b,a); zu[].zero, zu[].one (see cdisort.h)
       zz     :  Beam source vectors in eq. SS(19), STWL(24b)
       zzg    :  Beam source vectors in eq. KS(10)for a general source constant over a layer
       plk    :  Thermal source vectors z0,z1 by solving eq. SS(16),
                 Y-sub-zero,Y-sub-one in STWL(26)
       zbeam  :  Incident-beam source vectors


    O U T P U T    V A R I A B L E S:

       uum    :  Azimuthal components of the intensity in eq. STWJ(5),
                 STWL(6)

    I N T E R N A L    V A R I A B L E S:

       bnddir :  Direct intensity down at the bottom boundary
       bnddfu :  Diffuse intensity down at the bottom boundary
       bndint :  Intensity attenuated at both boundaries, STWJ(25-6)
       dtau   :  Optical depth of a computational layer
       lyrend :  End layer of integration
       lyrstr :  Start layer of integration
       palint :  Intensity component from parallel beam
       plkint :  Intensity component from planck source
       wk     :  Scratch vector for saving exp evaluations

       All the exponential factors (exp1, expn,... etc.)
       come from the substitution of constants of integration in
       eq. SC(12) into eqs. S1(8-9).  They all have negative
       arguments so there should never be overflow problems.

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_user_intensities(disort_state   *ds,
                        double          bplanck,
                        double         *cmu,
                        double         *cwt,
                        double          delm0,
                        double         *dtaucpr,
                        double         *emu,
                        double         *expbea,
                        double         *gc,
                        double         *gu,
                        double         *kk,
                        int            *layru,
                        double         *ll,
                        int             lyrcut,
                        int             mazim,
                        int             ncut,
                        int             nn,
                        double         *rmu,
                        double         *taucpr,
                        double          tplanck,
                        double         *utaupr,
                        double         *wk,
			disort_triplet *zbu,
                        double         *zbeam,
			disort_pair    *zbeamsp,
                        double         *zbeama,
                        double         *zgu,
                        disort_pair    *zu,
                        double         *zz,
                        double         *zzg,
                        disort_pair    *plk,
                        double         *uum)
{
  register int
    negumu,
    iq,iu,jq,lc,lu,lyrend,lyrstr,lyu;
  double
    alfa,bnddfu,bnddir,bndint,
    denom,dfuint,dtau,dtau1,dtau2,
    exp0=0,exp1=0,exp2=0,expn,
    f0n,f1n,fact,genint,
    palint,plkint,sgn;

  /*
   * Incorporate constants of integration into interpolated eigenvectors
   */
  for (lc = 1; lc <= ncut; lc++) {
    for (iq = 1; iq <= ds->nstr; iq++) {
      for (iu = 1; iu <= ds->numu; iu++) {
        GU(iu,iq,lc) *= LL(iq,lc);
      }
    }
  }

  /*
   * Loop over levels at which intensities are desired ('user output levels')
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    if (ds->bc.fbeam > 0.) {
      exp0 = exp(-UTAUPR(lu)/ds->bc.umu0);
    }
    lyu = LAYRU(lu);
    /*
     * Loop over polar angles at which intensities are desired
     */
    for (iu = 1; iu <= ds->numu; iu++) {
      if (lyrcut && lyu > ncut) {
        continue;
      }
      negumu = (UMU(iu) < 0.);
      if (negumu) {
        lyrstr = 1;
        lyrend = lyu-1;
        sgn    = -1.;
      }
      else {
        lyrstr = lyu+1;
        lyrend = ncut;
        sgn    = 1.;
      }

      /*
       * For downward intensity, integrate from top to LYU-1 in eq. S1(8); for upward,
       * integrate from bottom to LYU+1 in S1(9)
       */
      genint = 0.;
      palint = 0.;
      plkint = 0.;
      for (lc = lyrstr; lc <= lyrend; lc++) {
        dtau = DTAUCPR(lc);
        exp1 = exp((UTAUPR(lu)-TAUCPR(lc-1))/UMU(iu));
        exp2 = exp((UTAUPR(lu)-TAUCPR(lc  ))/UMU(iu));

        if (ds->flag.planck && mazim == 0) {
          /*
           * Eqs. STWL(36b,c, 37b,c)
           */
          f0n     = sgn*(exp1-exp2);
          f1n     = sgn*((TAUCPR(lc-1)+UMU(iu))*exp1
                        -(TAUCPR(lc  )+UMU(iu))*exp2);
          plkint += Z0U(iu,lc)*f0n+Z1U(iu,lc)*f1n;
        }

        if (ds->bc.fbeam > 0.) {
	  if ( ds->flag.spher == TRUE ) {
	    denom  =  sgn*1.0/(ZBAU(iu,lc)*UMU(iu)+1.0);
	    palint += (ZB0U(iu,lc)*denom*(exp(-ZBAU(iu,lc)*TAUCPR(lc-1)) *exp1
					  -exp(-ZBAU(iu,lc)*TAUCPR(lc)) *exp2 )
		       +ZB1U(iu,lc)*denom*((TAUCPR(lc-1)+sgn*denom*UMU(iu))
					   *exp(-ZBAU(iu,lc)*TAUCPR(lc-1)) *exp1
					   -(TAUCPR(lc)+sgn*denom*UMU(iu) )
					   *exp(-ZBAU(iu,lc)*TAUCPR(lc))*exp2));
	  }
	  else {
	    denom = 1.+UMU(iu)/ds->bc.umu0;
	    if (fabs(denom) < 0.0001) {
	      /*
	       * L'Hospital limit
	       */
	      expn = (dtau/ds->bc.umu0)*exp0;
	    }
	    else {
	      expn = (exp1*EXPBEA(lc-1)
		      -exp2*EXPBEA(lc  ))*sgn/denom;
	    }
	    palint += ZBEAM(iu,lc)*expn;
	  }
        }
	if ( ds->flag.general_source ) {
          genint += ZGU(iu,lc)*sgn*(exp1-exp2);
	}
        /*
         * KK is negative
         */
        for (iq = 1; iq <= nn; iq++) {
          WK(iq) = exp(KK(iq,lc)*dtau);
          denom  = 1.+UMU(iu)*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            /*
             * L'Hospital limit
             */
            expn = (dtau/UMU(iu))*exp2;
          }
          else {
            expn = sgn*(exp1*WK(iq)-exp2)/denom;
          }
          palint += GU(iu,iq,lc)*expn;
        }

        /*
         * KK is positive
         */
        for (iq = nn+1; iq <= ds->nstr; iq++) {
          denom = 1.+UMU(iu)*KK(iq,lc);
          if (fabs(denom) < 0.0001) {
            /*
             * L'Hospital limit
             */
            expn = -(dtau/UMU(iu))*exp1;
          }
          else {
            expn = sgn*(exp1-exp2*WK(ds->nstr+1-iq))/denom;
          }
          palint += GU(iu,iq,lc)*expn;
        }
      }

      /*
       * Calculate contribution from user output level to next computational level
       */
      dtau1 = UTAUPR(lu)-TAUCPR(lyu-1);
      dtau2 = UTAUPR(lu)-TAUCPR(lyu  );

      if ((fabs(dtau1) >= 1.e-6 || !negumu) && (fabs(dtau2) >= 1.e-6 ||  negumu)) {
        if(negumu) {
          exp1 = exp(dtau1/UMU(iu));
        }
        else {
          exp2 = exp(dtau2/UMU(iu));
        }
        if (ds->bc.fbeam > 0.) {
	  if ( ds->flag.spher == TRUE ) {
	    if ( negumu ) {
	      expn = exp1;
	      alfa = ZBAU(iu,lyu);
	      denom = (-1.0/(alfa*UMU(iu)+1.));
	      palint += ZB0U(iu,lyu)*denom*(-exp(-alfa*UTAUPR(lu))
					    + expn*exp(-alfa*TAUCPR(lyu-1)))
		+ZB1U(iu,lyu)*denom*( -(UTAUPR(lu)-UMU(iu)*denom)*exp(-alfa*UTAUPR(lu))
				      +(TAUCPR(lyu-1)-UMU(iu)*denom)*expn*exp(-alfa*TAUCPR(lyu-1)));
	    }
	    else {
	      expn = exp2;
	      alfa = ZBAU(iu,lyu);
	      denom = (1.0/(alfa*UMU(iu)+1.0));
	      palint += ZB0U(iu,lyu)*denom*(exp(-alfa*UTAUPR(lu))
					    -exp(-alfa*TAUCPR(lyu))*expn)
		+ZB1U(iu,lyu)*denom*( (UTAUPR(lu) +UMU(iu)*denom)*exp(-alfa*UTAUPR(lu))
				      -(TAUCPR(lyu)+UMU(iu)*denom)*exp(-alfa*TAUCPR(lyu))*expn );
	    }
	  }
	  else {
	    denom = 1.+UMU(iu)/ds->bc.umu0;
	    if (fabs(denom) < 0.0001) {
	      expn = (dtau1/ds->bc.umu0)*exp0;
	    }
	    else if (negumu) {
	      expn = (exp0-EXPBEA(lyu-1)*exp1)/denom;
	    }
	    else {
	      expn = (exp0-EXPBEA(lyu  )*exp2)/denom;
	    }
	    palint += ZBEAM(iu,lyu)*expn;
	  }
        }
	if ( ds->flag.general_source ) {
          if (negumu) {
            expn = exp1;
          }
          else {
            expn = exp2;
          }
          genint += ZGU(iu,lyu)*(1.-expn);
	}
        /*
         * KK is negative
         */
        dtau = DTAUCPR(lyu);
        for (iq = 1; iq <= nn; iq++) {
          denom = 1.+UMU(iu)*KK(iq,lyu);
          if (fabs(denom) < 0.0001) {
            expn = -dtau2/UMU(iu)*exp2;
          }
          else if (negumu) {
            expn = (exp(-KK(iq,lyu)*dtau2)
                   -exp( KK(iq,lyu)*dtau )*exp1)/denom;
          }
          else {
            expn = (exp(-KK(iq,lyu)*dtau2)-exp2)/denom;
          }
          palint += GU(iu,iq,lyu)*expn;
        }

        /*
         * KK is positive
         */
        for (iq = nn+1; iq <= ds->nstr; iq++) {
          denom = 1.+UMU(iu)*KK(iq,lyu);
          if (fabs(denom) < 0.0001) {
            expn = -(dtau1/UMU(iu))*exp1;
          }
          else if (negumu) {
            expn = (exp(-KK(iq,lyu)*dtau1)-exp1)/denom;
          }
          else {
            expn = (exp(-KK(iq,lyu)*dtau1)
                   -exp(-KK(iq,lyu)*dtau )*exp2)/denom;
          }
          palint += GU(iu,iq,lyu)*expn;
        }

        if (ds->flag.planck && mazim == 0) {
          /*
           * Eqs. STWL (35-37) with tau-sub-n-1 replaced by tau for upward, and
           * tau-sub-n replaced by tau for downward directions
           */
          if (negumu) {
            expn = exp1;
            fact = TAUCPR(lyu-1)+UMU(iu);
          }
          else {
            expn = exp2;
            fact = TAUCPR(lyu  )+UMU(iu);
          }
          f0n     = 1.-expn;
          f1n     = UTAUPR(lu)+UMU(iu)-fact*expn;
          plkint += Z0U(iu,lyu)*f0n+Z1U(iu,lyu)*f1n;
        }
      }

      /*
       * Calculate intensity components attenuated at both boundaries.
       * NOTE: no azimuthal intensity component for isotropic surface
       */
      bndint = 0.;
      if (negumu && mazim == 0) {
        bndint = (ds->bc.fisot+tplanck)*exp(UTAUPR(lu)/UMU(iu));
      }
      else if (!negumu) {
        if (lyrcut || ( ds->flag.lamber && mazim > 0 ) ) {
          UUM(iu,lu) = palint+plkint;
          continue;
        }

        for (jq = nn+1; jq <= ds->nstr; jq++) {
          WK(jq) = exp(-KK(jq,ds->nlyr)*DTAUCPR(ds->nlyr));
        }
        bnddfu = 0.;
        for (iq = nn; iq >= 1; iq--) {
          dfuint = 0.;
          for (jq = 1; jq <= nn; jq++) {
            dfuint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr);
          }
          for (jq= nn+1; jq <= ds->nstr; jq++) {
            dfuint += GC(iq,jq,ds->nlyr)*LL(jq,ds->nlyr)*WK(jq);
          }
          if (ds->bc.fbeam > 0.) {
	    if ( ds->flag.spher == TRUE ) {
	      dfuint += exp(-ZBEAMA(ds->nlyr)*TAUCPR(ds->nlyr)) *
		(ZBEAM0(iq,ds->nlyr)+ZBEAM1(iq,ds->nlyr)*TAUCPR(ds->nlyr));
	    }
	    else {
	      dfuint += ZZ(iq,ds->nlyr)*EXPBEA(ds->nlyr);
	    }
          }
	  if ( ds->flag.general_source ) {
	    dfuint += ZZG(iq,ds->nlyr);
	  }
          dfuint += delm0*(ZPLK0(iq,ds->nlyr)+ZPLK1(iq,ds->nlyr)*TAUCPR(ds->nlyr));
          bnddfu += (1.+delm0)*RMU(iu,nn+1-iq)*CMU(nn+1-iq)*CWT(nn+1-iq)*dfuint;
        }
        bnddir = 0.;
        if (ds->bc.fbeam > 0. || ds->bc.umu0 >0.) {
          bnddir = ds->bc.umu0*ds->bc.fbeam/M_PI*RMU(iu,0)*EXPBEA(ds->nlyr);
        }
        bndint = (bnddfu+bnddir+delm0*EMU(iu)*bplanck+ds->bc.fluor)*exp((UTAUPR(lu)-TAUCPR(ds->nlyr))/UMU(iu));
      }
      UUM(iu,lu) = palint+plkint+bndint+genint;
    }
  }

  return;
}

/*============================= end of c_user_intensities() =============*/
