// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_disort_set() ==========================*/

/*
    Perform miscellaneous setting-up operations

    I N P U T  V A R I A B L E S

       ds         Disort state variables
       deltam
       tauc

    O U T P U T     V A R I A B L E S:

       If ds->flag.usrtau is FALSE
       ds->ntau
       ds->utau

       If ds->flag.usrang is FALSE
       ds->numu
       ds->umu

       cmu,cwt     computational polar angles and corresponding quadrature weights
       dtaucpr
       expbea      transmission of direct beam
       flyr        separated fraction in delta-m method
       gl          phase function legendre coefficients multiplied by (2l+1) and single-scatter albedo
       layru       computational layer in which utau falls
       lyrcut      flag as to whether radiation will be zeroed below layer ncut
       ncut        computational layer where absorption optical depth first exceeds  abscut
       nn          ds->nstr/2
       oprim       delta-m-scaled single-scatter albedo
       taucpr      delta-m-scaled optical depth
       utaupr      delta-m-scaled version of  utau

   Called by- c_disort
   Calls- c_gaussian_quadrature, c_errmsg

 ---------------------------------------------------------------------*/

void c_disort_set(disort_state *ds,
                  double       *ch,
                  double       *chtau,
                  double       *cmu,
                  double       *cwt,
                  int           deltam,
                  double       *dtaucpr,
                  double       *expbea,
                  double       *flyr,
                  double       *gl,
                  int          *layru,
                  int          *lyrcut,
                  int          *ncut,
                  int          *nn,
                  int          *corint,
                  double       *oprim,
                  double       *tauc,
                  double       *taucpr,
                  double       *utaupr,
                  emission_func_t emi_func)
{
  register int
    iq,iu,k,lc,lu;
  const double
    abscut = 10.;
  double
    abstau,chtau_tmp,f,taup,zenang;

  if (!ds->flag.usrtau) {
   /*
    * Set output levels at computational layer boundaries
    */
    for (lc = 0;  lc <= ds->ntau-1; lc++) {
      UTAU(lc+1) = TAUC(lc);
    }
  }

  /*
   * Apply delta-M scaling and move description of computational layers to local variables
   */
  TAUCPR(0) = 0.;
  abstau    = 0.;
  for (lc = 1; lc <= ds->nlyr; lc++) {
    PMOM(0,lc)  = 1.;
    if (abstau < abscut) {
      *ncut = lc;
    }
    abstau += (1.-SSALB(lc))*DTAUC(lc);
    if (!deltam) {
      OPRIM(lc)   = SSALB(lc);
      DTAUCPR(lc) = DTAUC(lc);
      TAUCPR(lc)  = TAUC(lc);
      for (k = 0; k <= ds->nstr-1; k++) {
        GL(k,lc)  = (double)(2*k+1)*OPRIM(lc)*PMOM(k,lc);
      }
      f = 0.;
    }
    else {
      /*
       * Do delta-M transformation
       */
      f           = PMOM(ds->nstr,lc);
      OPRIM(lc)   = SSALB(lc)*(1.-f)/(1.-f*SSALB(lc));
      DTAUCPR(lc) = (1.-f*SSALB(lc))*DTAUC(lc);
      TAUCPR(lc)  = TAUCPR(lc-1)+DTAUCPR(lc);
      for (k = 0; k <= ds->nstr-1; k++) {
        GL(k,lc)  = (double)(2*k+1)*OPRIM(lc)*(PMOM(k,lc)-f)/(1.-f);
      }
    }

    FLYR(lc)   = f;
  }

  /*
   * Calculate Chapman function if spherical geometry, set expbea and
   * ch for beam source.
   */
  if( (ds->flag.ibcnd == GENERAL_BC && ds->bc.fbeam > 0.) ||
      (ds->flag.ibcnd == GENERAL_BC && ds->flag.general_source )) {

    CHTAU(0)  = 0.;
    EXPBEA(0) = 1.;
    zenang    = acos(ds->bc.umu0)/DEG;

    if( ds->flag.spher == TRUE && ds->bc.umu0 < 0. ) {
      EXPBEA(0) = exp(-c_chapman(1,0.,tauc,ds->nlyr,ds->zd,
				 ds->dtauc,zenang,ds->radius));
    }
    if ( ds->flag.spher == TRUE ) {
      for (lc = 1; lc <= *ncut; lc++) {
        taup        = TAUCPR(lc-1) + DTAUCPR(lc)/2.;
	/* Need Chapman function at top (0.0) and middle (0.5) of layer */
        CHTAU(lc  ) = c_chapman(lc, 0.,   taucpr,ds->nlyr,ds->zd,
				dtaucpr,zenang,ds->radius);
        chtau_tmp   = c_chapman(lc, 0.5,  taucpr,ds->nlyr,ds->zd,
				dtaucpr,zenang,ds->radius);
        CH(lc)      = taup/chtau_tmp;
        EXPBEA(lc)  = exp(-CHTAU(lc));
      }
    }
    else {
      for (lc = 1; lc <= *ncut; lc++) {
        CH(lc)     = ds->bc.umu0;
        EXPBEA(lc) = exp(-TAUCPR(lc)/ds->bc.umu0);
      }
    }
  }
  else {
    for (lc = 1; lc <= *ncut; lc++) {
      EXPBEA(lc) = 0.;
    }
  }

  /*
   * If no thermal emission, cut off medium below absorption optical depth = abscut ( note that
   * delta-M transformation leaves absorption optical depth invariant ).  Not worth the
   * trouble for one-layer problems, though.
   */
  *lyrcut = FALSE;
  if (abstau >= abscut && !ds->flag.planck && ds->flag.ibcnd != SPECIAL_BC && ds->nlyr > 1) {
    *lyrcut = TRUE;
  }
  if(!*lyrcut) *ncut = ds->nlyr;

  /*
   * Set arrays defining location of user output levels within delta-M-scaled computational mesh
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    for (lc = 1; lc < ds->nlyr; lc++) {
      if (UTAU(lu) >= TAUC(lc-1) && UTAU(lu) <= TAUC(lc)) {
        break;
      }
    }

    UTAUPR(lu) = UTAU(lu);
    if (deltam) {
      UTAUPR(lu) = TAUCPR(lc-1)+(1.-SSALB(lc)*FLYR(lc))*(UTAU(lu)-TAUC(lc-1));
    }
    LAYRU(lu) = lc;
  }

  /*
   * Calculate computational polar angle cosines and associated quadrature weights for Gaussian
   * quadrature on the interval (0,1) (upward)
   */
  *nn = ds->nstr/2;
  c_gaussian_quadrature(*nn,cmu,cwt);

  /*
   * Downward (neg) angles and weights
   */
  for (iq = 1; iq <= *nn; iq++) {
    CMU(iq+*nn) = -CMU(iq);
    CWT(iq+*nn) =  CWT(iq);
  }

  if (ds->flag.ibcnd == GENERAL_BC && ds->bc.fbeam > 0.) {
    /*
     * Compare beam angle to comput. angles
     */
    for (iq = 1; iq <= *nn; iq++) {
      if (fabs(ds->bc.umu0-CMU(iq))/fabs(ds->bc.umu0) < 1.e-4) {
        // suppress error msg by adding a small difference
        ds->bc.umu0 = (1. + 1.E-4)*CMU(iq);
        // c_errmsg("cdisort_set--beam angle=computational angle; change ds.nstr",DS_ERROR);
      }
    }
  }

  if (!ds->flag.usrang || ds->flag.onlyfl) {
    /*
     * Set output polar angles to computational polar angles
     */
    for (iu = 1; iu <= *nn; iu++) {
      UMU(iu) = -CMU(*nn+1-iu);
    }
    for (iu = *nn+1; iu <=ds->nstr; iu++) {
      UMU(iu) =  CMU(iu-*nn);
    }
  }

  if (ds->flag.usrang && ds->flag.ibcnd == SPECIAL_BC) {
    /*
     * Shift positive user angle cosines to upper locations and put negatives in lower locations
     */
    for (iu = 1; iu <= ds->numu/2; iu++) {
      UMU(iu+ds->numu/2) = UMU(iu);
    }
    for (iu = 1; iu <= ds->numu/2; iu++) {
      UMU(iu) = -UMU((ds->numu/2)+1-iu);
    }
  }

  return;
}

/*============================= end of c_disort_set() ===================*/
