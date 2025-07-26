// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_solve0() ==============================*/

/*
        Construct right-hand side vector B for general boundary
        conditions STWJ(17) and solve system of equations obtained
        from the boundary conditions and the continuity-of-
        intensity-at-layer-interface equations.
        Thermal emission contributes only in azimuthal independence.

    I N P U T      V A R I A B L E S:

       ds       :  Disort input variables
       bdr      :  Surface bidirectional reflectivity
       bem      :  Surface bidirectional emissivity
       bplanck  :  Bottom boundary thermal emission
       cband    :  Left-hand side matrix of linear system eq. SC(5),
                   scaled by eq. SC(12); in banded form required
                   by LINPACK solution routines
       cmu,cwt  :  Abscissae, weights for Gauss quadrature
                   over angle cosine
       expbea   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
       lyrcut   :  Logical flag for truncation of computational layers
       mazim    :  Order of azimuthal component
       ncol     :  Number of columns in CBAND
       nn       :  Order of double-Gauss quadrature (NSTR/2)
       ncut     :  Total number of computational layers considered
       tplanck  :  Top boundary thermal emission
       taucpr   :  Cumulative optical depth (delta-M-scaled)
       zz       :  Beam source vectors in eq. SS(19), STWL(24b)
       zzg      :  Beam source vectors in eq. KS(10)for a general source constant over a layer
       plk      :  Thermal source vectors z0,z1 by solving eq. SS(16), Y0,Y1 in STWL(26b,a);
                   plk[].zero, plk[].one (see cdisort.h)

    O U T P U T     V A R I A B L E S:

       b        :  Right-hand side vector of eq. SC(5) going into
                   sgbsl; returns as solution vector of eq. SC(12),
                   constants of integration without exponential term
      ll        :  Permanent storage for B, but re-ordered

   I N T E R N A L    V A R I A B L E S:

       ipvt     :  Integer vector of pivot indices
       it       :  Pointer for position in  B
       ncd      :  Number of diagonals below or above main diagonal
       rcond    :  Indicator of singularity for cband
       z        :  Scratch array required by sgbco

   Called by- c_disort
   Calls- c_sgbco, c_errmsg, c_sgbsl
 +-------------------------------------------------------------------*/

void c_solve0(disort_state *ds,
              double       *b,
              double       *bdr,
              double       *bem,
              double        bplanck,
              double       *cband,
              double       *cmu,
              double       *cwt,
              double       *expbea,
              int          *ipvt,
              double       *ll,
              int           lyrcut,
              int           mazim,
              int           ncol,
              int           ncut,
              int           nn,
              double        tplanck,
              double       *taucpr,
              double       *z,
              disort_pair  *zbeamsp,
	      double       *zbeama,
              double       *zz,
              double       *zzg,
              disort_pair  *plk)
{
  register int
    ipnt,iq,it,jq,lc,ncd;
  double
    rcond,sum,diff;

  memset(b,0,ds->nstr*ds->nlyr*sizeof(double));

  /*
   * Construct B, STWJ(20a,c) for parallel beam+bottom
   * reflection+thermal emission at top and/or bottom
   */
  if (mazim > 0 && ( ds->bc.fbeam > 0.  || ds->flag.general_source) ) {
    /*
     * Azimuth-dependent case (never called if FBEAM = 0)
     */
    if ( lyrcut == TRUE || ds->flag.lamber == TRUE ) {
      /*
       * No azimuthal-dependent intensity for Lambert surface; no
       * intensity component for truncated bottom layer
       */
      for (iq = 1; iq <= nn; iq++) {
        /*
         * Top boundary
         */
	if ( ds->flag.spher == TRUE ) {
	  B(iq) = - ZBEAM0(nn+1-iq,1);
	}
	else {
	  B(iq) = - ZZ(nn+1-iq,1);
	}
	if ( ds->flag.general_source == TRUE ) {
	  B(iq) -= ZZG(nn+1-iq,1);
	  //aky	  B(iq) = B(iq) - ZZG(nn+1-iq,1);
	}
        /*
         * Bottom boundary
         */
	if ( ds->flag.spher == TRUE ) {
	  B(ncol-nn+iq) = - exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
	    (ZBEAM0(iq+nn,ncut) + ZBEAM1(iq+nn,ncut)*TAUCPR(ncut));
	}
	else {
	  B(ncol-nn+iq) = - ZZ(iq+nn,ncut)*EXPBEA(ncut);
	}
	if ( ds->flag.general_source == TRUE ) {
	  B(ncol-nn+iq) -=  ZZG(iq+nn,ncut);
	  //aky	  B(ncol-nn+iq) = B(ncol-nn+iq)  - ZZG(iq+nn,ncut);
	}
      }
    }
    else {
      for (iq = 1; iq <= nn; iq++) {
	if ( ds->flag.spher == TRUE ) {
	  B(iq) = - ZBEAM0(nn+1-iq,1);
	}
	else {
	  B(iq) = - ZZ(nn+1-iq,1);
	}
	if ( ds->flag.general_source == TRUE ) {
	  B(iq) -= ZZG(nn+1-iq,1);
	  //aky	  B(iq) = B(iq) - ZZG(nn+1-iq,1);
	}
	if ( ds->flag.spher == TRUE ) {
	  c_errmsg("solve0--BDR not implemented for pseudo-spherical geometry",
		   DS_WARNING);
	}
	else {
	  sum   = 0.;
	  for (jq = 1; jq <= nn; jq++) {
	    sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*ZZ(nn+1-jq,ncut)*EXPBEA(ncut);
	  }
	  B(ncol-nn+iq) = sum;
	  if ( ds->flag.general_source == TRUE ) {
	    sum   = 0.;
	    for (jq = 1; jq <= nn; jq++) {
	      sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*ZZG(nn+1-jq,ncut);
	    }
	    B(ncol-nn+iq) += sum;
	  }
	}
        if (ds->bc.fbeam > 0.) {
	  if ( ds->flag.spher == TRUE ) {
	    c_errmsg("solve0--BDR not implemented for pseudo-spherical geometry",
		     DS_WARNING)  ;
	  }
	  else {
	    B(ncol-nn+iq) += (BDR(iq,0)*ds->bc.umu0*ds->bc.fbeam/
			      M_PI-ZZ(iq+nn,ncut))*EXPBEA(ncut);
	  }
        }
	if ( ds->flag.general_source == TRUE ) {
	    B(ncol-nn+iq) += -ZZG(iq+nn,ncut);
	}
      }
    }
    /*
     * Continuity condition for layer interfaces of eq. STWJ(20b)
     */
    it = nn;
    diff = 0;
    for (lc = 1; lc <= ncut-1; lc++) {
      for (iq = 1; iq <= ds->nstr; iq++) {
	if ( ds->flag.general_source == TRUE ) {
	  diff = (ZZG(iq,lc+1)-ZZG(iq,lc));
	}
	if ( ds->flag.spher == TRUE ) {
	  B(++it) = exp(-ZBEAMA(lc+1)*TAUCPR(lc))*
	    (ZBEAM0(iq,lc+1)+ZBEAM1(iq,lc+1)*TAUCPR(lc))
	    -  exp(-ZBEAMA(lc)*TAUCPR(lc))*
	    (ZBEAM0(iq,lc)+ZBEAM1(iq,lc)*TAUCPR(lc))
	    + diff;
	}
	else {
	  B(++it) = (ZZ(iq,lc+1)-ZZ(iq,lc))*EXPBEA(lc)  + diff;
	}
      }
    }
  }
  else {
    /*
     * Azimuth-independent case
     */
    if (ds->bc.fbeam == 0. && ds->flag.general_source == FALSE ) {
      for (iq = 1; iq <= nn; iq++) {
        /*
         * Top boundary
         */
        B(iq) = -ZPLK0(nn+1-iq,1)+ds->bc.fisot+tplanck;
      }
      if ( lyrcut == TRUE ) {
        /*
         * No intensity component for truncated bottom layer
         */
        for (iq = 1; iq <= nn; iq++) {
          /*
           * Bottom boundary
           */
          B(ncol-nn+iq) = -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
        }
      }
      else {
        for (iq = 1; iq <= nn; iq++) {
          sum = 0.;
          for (jq = 1; jq <= nn; jq++) {
            sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*
	      (ZPLK0(nn+1-jq,ncut)+ZPLK1(nn+1-jq,ncut)*TAUCPR(ncut));
          }
          B(ncol-nn+iq) = 2.*sum+BEM(iq)*bplanck-
	    ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
        }
      }
      /*
       * Continuity condition for layer interfaces, STWJ(20b)
       */
      it = nn;
      for (lc = 1; lc <= ncut-1; lc++) {
        for (iq = 1; iq <= ds->nstr; iq++) {
          B(++it) = ZPLK0(iq,lc+1)-ZPLK0(iq,lc)+
	    (ZPLK1(iq,lc+1)-ZPLK1(iq,lc))*TAUCPR(lc);
        }
      }
    }
    else {
      if ( ds->flag.spher == TRUE ) {
	for (iq = 1; iq <= nn; iq++)
	  B(iq) = -ZBEAM0(nn+1-iq,1)-ZPLK0(nn+1-iq,1)+ds->bc.fisot+tplanck;
      }
      else {
	for (iq = 1; iq <= nn; iq++)
	  B(iq) = -ZZ(nn+1-iq,1)-ZPLK0(nn+1-iq,1)+ds->bc.fisot+tplanck;
      }
      if ( ds->flag.general_source == TRUE ) {
	for (iq = 1; iq <= nn; iq++)
	  B(iq) -= ZZG(nn+1-iq,1);
	//aky	  B(iq) = B(iq) - ZZG(nn+1-iq,1);
      }
      if (lyrcut) {
	if ( ds->flag.spher == TRUE ) {
	  for (iq = 1; iq <= nn; iq++) {
	    B(ncol-nn+iq) = -exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
	      (ZBEAM0(iq+nn,ncut)+ ZBEAM1(iq+nn,ncut)*TAUCPR(ncut))
	      -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
	  }
	}
	else {
	  for (iq = 1; iq <= nn; iq++) {
	    B(ncol-nn+iq) = -ZZ(iq+nn,ncut)*EXPBEA(ncut)
	      -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut);
	  }
	}
	if ( ds->flag.general_source == TRUE ) {
	  for (iq = 1; iq <= nn; iq++)
	    B(ncol-nn+iq) -= ZZG(iq+nn,ncut);
	  //aky	    B(ncol-nn+iq) = B(ncol-nn+iq) - ZZG(iq+nn,ncut);
	}
      }
      else {
	if ( ds->flag.spher == TRUE ) {
	  for (iq = 1; iq <= nn; iq++) {
	    sum = 0.;
	    for (jq = 1; jq <= nn; jq++) {
	      sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*
		( exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
		  (ZBEAM0(nn+1-jq,ncut)+ZBEAM1(nn+1-jq,ncut)*TAUCPR(ncut))
		  + ZZG(nn+1-jq,ncut)
		  + ZPLK0(nn+1-jq,ncut)+ZPLK1(nn+1-jq,ncut)*TAUCPR(ncut));
	    }
	    B(ncol-nn+iq) = 2.0*sum +
	      ( BDR(iq,0)*ds->bc.umu0*ds->bc.fbeam/M_PI) *EXPBEA(ncut)
	      -  exp(-ZBEAMA(ncut)*TAUCPR(ncut))*
	      (ZBEAM0(iq+nn,ncut)+ZBEAM1(iq+nn,ncut)*TAUCPR(ncut))
	      - ZZG(iq+nn,ncut)
	      + BEM(iq)*bplanck
	      -ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut)
	      +ds->bc.fluor;
	  }
	}
	else {
	  for (iq = 1; iq <= nn; iq++) {
	    sum = 0.;
	    for (jq = 1; jq <= nn; jq++) {
	      sum += CWT(jq)*CMU(jq)*BDR(iq,jq)*
		(ZZ(nn+1-jq,ncut)*EXPBEA(ncut)+ZPLK0(nn+1-jq,ncut)
		 + ZZG(nn+1-jq,ncut)
		 +ZPLK1(nn+1-jq,ncut)*TAUCPR(ncut));
	    }
	    B(ncol-nn+iq) = 2.*sum+
	      (BDR(iq,0)*ds->bc.umu0*ds->bc.fbeam/M_PI-ZZ(iq+nn,ncut))
	      *EXPBEA(ncut)
	      - ZZG(iq+nn,ncut)
	      +BEM(iq)*bplanck-ZPLK0(iq+nn,ncut)-ZPLK1(iq+nn,ncut)*TAUCPR(ncut)
	      +ds->bc.fluor;
	  }
	}
      }
      it = nn;
      if ( ds->flag.spher == TRUE ) {
	for (lc = 1; lc <= ncut-1; lc++) {
	  for (iq = 1; iq <= ds->nstr; iq++) {
	    B(++it) = exp(-ZBEAMA(lc+1)*TAUCPR(lc))*
	      (ZBEAM0(iq,lc+1)+ZBEAM1(iq,lc+1)*TAUCPR(lc))
	      -exp(-ZBEAMA(lc)*TAUCPR(lc))*
	      (ZBEAM0(iq,lc)+ZBEAM1(iq,lc)*TAUCPR(lc))
	      +ZZG(iq,lc+1)-ZZG(iq,lc)
	      +ZPLK0(iq,lc+1)-ZPLK0(iq,lc)+
	      (ZPLK1(iq,lc+1)-ZPLK1(iq,lc))*TAUCPR(lc);
	  }
	}
      }
      else {
	for (lc = 1; lc <= ncut-1; lc++) {
	  for (iq = 1; iq <= ds->nstr; iq++) {
	    B(++it) = (ZZ(iq,lc+1)-ZZ(iq,lc))*EXPBEA(lc)
	      +ZZG(iq,lc+1)-ZZG(iq,lc)
	      +ZPLK0(iq,lc+1)-ZPLK0(iq,lc)
	      +(ZPLK1(iq,lc+1)-ZPLK1(iq,lc))*TAUCPR(lc);
	  }
	}
      }
    }
  }

  /*
   * Find L-U (lower/upper triangular) decomposition of band matrix
   * CBAND and test if it is nearly singular (note: CBAND is
   * destroyed) (CBAND is in LINPACK packed format)
   */
  rcond = 0.;
  ncd   = 3*nn-1;
  c_sgbco(cband,(9*(ds->nstr/2)-2),ncol,ncd,ncd,ipvt,&rcond,z);

  if (1.+rcond == 1.) {
    c_errmsg("solve0--sgbco says matrix near singular",DS_WARNING);
  }

  /*
   * Solve linear system with coeff matrix CBAND and R.H. side(s) B
   * after CBAND has been L-U decomposed. Solution is returned in B.
   */

  c_sgbsl(cband,(9*(ds->nstr/2)-2),ncol,ncd,ncd,ipvt,b,0);

  /*
   * Zero CBAND (it may contain 'foreign' elements upon returning from
   * LINPACK); necessary to prevent errors
   */
  memset(cband,0,(9*(ds->nstr/2)-2)*(ds->nstr*ds->nlyr)*sizeof(double));

  for (lc = 1; lc <= ncut; lc++) {
    ipnt = lc*ds->nstr-nn;
    for (iq = 1; iq <= nn; iq++) {
      LL(nn-iq+1,lc) = B(ipnt-iq+1);
      LL(nn+iq,  lc) = B(ipnt+iq  );
    }
  }

  return;
}
/*============================= end of c_solve0() =======================*/
