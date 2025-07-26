// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_upbeam_pseudo_spherical() =============*/

/*

       Finds the particular solution of beam source KS(10-11)

     Routines called:  sgeco, sgesl

   I N P U T     V A R I A B L E S:

       cc     :  capital-c-sub-ij in Eq. SS(5)
       cmu    :  abscissae for gauss quadrature over angle cosine
       xb0    :  EXPansion of beam source function Eq. KS(7)
       xb1    :  EXPansion of beam source function Eq. KS(7)
       xba    :  EXPansion of beam source function Eq. KS(7)
       (remainder are 'disort' input variables)

    O U T P U T    V A R I A B L E S:

       zbs0     :  solution vectors z-sub-zero of Eq. KS(10-11)
       zbs1     :  solution vectors z-sub-one  of Eq. KS(10-11)
       zbsa     :  alfa coefficient in Eq. KS(7)
       zbeam0,  :  permanent storage for -zbs0,zbs1,zbsa-, but rD-ordered
        zbeam1,
        zbeama

   I N T E R N A L    V A R I A B L E S:

       array  :  coefficient matrix in left-hand side of Eq. KS(10)
       ipvt   :  integer vector of pivot indices required by *linpack*
       wk     :  scratch array required by *linpack*

   Called by- c_disort
   Calls- c_sgeco, c_errmsg, c_sgesl
 -------------------------------------------------------------------*/

#undef  ARRAY
#define ARRAY(iq,jq) array[iq-1+(jq-1)*ds->nstr]

void c_upbeam_pseudo_spherical(disort_state *ds,
			       int           lc,
			       double       *array,
			       double       *cc,
			       double       *cmu,
			       int          *ipvt,
			       int           nn,
			       double       *wk,
			       disort_pair  *xb,
			       double       *xba,
			       disort_pair  *zbs,
			       double       *zbsa,
			       disort_pair  *zbeamsp,
			       double       *zbeama)
{

  register int
    iq,jq;
  double
    rcond,rmin;


  for (iq = 1; iq <= ds->nstr; iq++) {
    for (jq = 1; jq <= ds->nstr; jq++) {
      ARRAY(iq,jq) = -CC(iq,jq);
    }
    ARRAY(iq,iq) += 1.+XBA(lc)*CMU(iq);
    *zbsa     = XBA(lc);
    ZBS1(iq) = XB1(iq,lc);
  }

  /*
   * Find L-U (lower/upper triangular) decomposition of ARRAY and see
   * if it is nearly singular
   * (NOTE: ARRAY is altered)
   */

  rcond = 0.;
  c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);

  if (1.+rcond == 1.) {
    c_errmsg("upbeam_pseudo_spherical--sgeco says matrix near singular",
	     DS_WARNING);
  }

  rmin = 1.0e-4;
  if ( rcond < rmin ) {
    /*     Dither alpha if rcond to small   */
    if(XBA(lc) ==0.0)       XBA(lc)=0.000000005;

    XBA(lc) = XBA(lc) * 1.00000005;

    for (iq = 1; iq <= ds->nstr; iq++) {
      for (jq = 1; jq <= ds->nstr; jq++) {
	ARRAY(iq,jq) = -CC(iq,jq);
      }
      ARRAY(iq,iq) += 1.0+XBA(lc)*CMU(iq);
      *zbsa     = XBA(lc);
      ZBS1(iq) = XB1(iq,lc);
    }
    /*     Solve linear equations KS(10-11) with dithered alpha */
    rcond = 0.;
    c_sgeco(array,ds->nstr,ds->nstr,ipvt,&rcond,wk);
    if (1.+rcond == 1.) {
      c_errmsg("upbeam_pseudo_spherical--sgeco says matrix near singular",
	       DS_WARNING);
    }
  }

  for (iq = 1; iq <= ds->nstr; iq++)  WK(iq) = ZBS1(iq);
  c_sgesl( array, ds->nstr, ds->nstr, ipvt, wk, 0 );

  for (iq = 1; iq <= ds->nstr; iq++) {
    ZBS1(iq) = WK(iq);
    ZBS0(iq) = XB0(iq,lc) + CMU(iq) * ZBS1(iq);
  }

  for (iq = 1; iq <= ds->nstr; iq++)  WK(iq) = ZBS0(iq);
  c_sgesl( array, ds->nstr, ds->nstr, ipvt, wk, 0 );
  for (iq = 1; iq <= ds->nstr; iq++)  ZBS0(iq) = WK(iq);

  /*   ... and now some index gymnastic for the inventive ones...  */

  ZBEAMA(lc)            = *zbsa;
  for (iq = 1; iq <= nn; iq++) {
    ZBEAM0( iq+nn, lc )   = ZBS0( iq );
    ZBEAM1( iq+nn, lc )   = ZBS1( iq );
    ZBEAM0( nn+1-iq, lc ) = ZBS0( iq+nn );
    ZBEAM1( nn+1-iq,lc )  = ZBS1( iq+nn );
  }

 return;

}


/*============================= end of c_upbeam_pseudo_spherical() ======*/
