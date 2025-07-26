// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_set() ===========================*/

/*
 Perform miscellaneous setting-up operations

 Routines called: c_errmsg

 Input :  ds         'Disort' input variables

 Output:  ntau,utau  If ds->flag.usrtau = FALSE
          bplanck    Intensity emitted from bottom boundary
          ch         The Chapman factor
          cmu        Computational polar angle
          expbea     Transmission of direct beam
          flyr       Truncated fraction in delta-m method
          layru      Computational layer in which utau falls
          lyrcut     Flag as to whether radiation will be zeroed below layer ncut
          ncut       Computational layer where absorption optical depth first exceeds abscut
          nn         nstr/2 = 1
          nstr       No.of streams (=2)
          oprim      Delta-m-scaled single-scatter albedo
          pkag,c     Planck function in each layer
          taucpr     Delta-m-scaled optical depth
          tplanck    Intensity emitted from top boundary
          utaupr     Delta-m-scaled version of utau

 Internal Variables
          abscut     Absorption optical depth, medium is cut off below this depth
          tempc      Temperature at center of layer, assumed to be average of
                     layer boundary temperatures
  ---------------------------------------------------------------------*/

void c_twostr_set(disort_state *ds,
                  double       *bplanck,
                  double       *ch,
                  double       *chtau,
                  double       *cmu,
                  int           deltam,
                  double       *dtaucpr,
                  double       *expbea,
                  double       *flyr,
                  double       *gg,
                  double       *ggprim,
                  int          *layru,
                  int          *lyrcut,
                  int          *ncut,
                  int          *nn,
                  double       *oprim,
                  double       *pkag,
                  double       *pkagc,
                  double        radius,
                  double       *tauc,
                  double       *taucpr,
                  double       *tplanck,
                  double       *utaupr,
                  emission_func_t emi_func)
{
  static int
    firstpass = TRUE;
  register int
    lc,lu,lev;
  double
    zenang,abstau,chtau_tmp,f,tempc,taup,
    abscut = 10.;

  if (firstpass) {
    firstpass = FALSE;
    ds->nstr  = 2;
    *nn       = ds->nstr/2;
  }

  if (!ds->flag.usrtau) {
    /*
     * Set output levels at computational layer boundaries
     */
    ds->ntau = ds->nlyr+1;
    for (lc = 0; lc <= ds->ntau-1; lc++) {
      UTAU(lc+1) = TAUC(lc);
    }
  }
  /*
   * Apply delta-m scaling and move description of computational layers to local variables
   */

  /*
   * NOTE: If not using swappablecalloc() to dynamically allocate memory, then need to zero-out
   *       taucpr, expbea, flyr, oprim here.
   */

  abstau = 0.;
  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (abstau < abscut) {
      *ncut = lc;
    }
    abstau += (1.-SSALB(lc))*DTAUC(lc);
    if (!deltam) {
      OPRIM(lc)   = SSALB(lc);
      TAUCPR(lc)  = TAUC(lc);
      f           = 0.;
      GGPRIM(lc)  = GG(lc);
      DTAUCPR(lc) = DTAUC(lc);
    }
    else {
     /*
      * Do delta-m transformation eqs. WW(20a,20b,14)
      */
      f           = SQR(GG(lc));
      TAUCPR(lc)  = TAUCPR(lc-1)+(1.-f*SSALB(lc))*DTAUC(lc);
      OPRIM(lc)   = SSALB(lc)*(1.-f)/(1.-f*SSALB(lc));
      GGPRIM(lc)  = (GG(lc)-f)/(1.-f);
      DTAUCPR(lc) = TAUCPR(lc)-TAUCPR(lc-1);
    }
    FLYR(lc) = f;
  }
  /*
   * If no thermal emission, cut off medium below absorption optical
   * depth = abscut (note that delta-m transformation leaves absorption
   * optical depth invariant). Not worth the trouble for one-layer problems, though.
   */
  *lyrcut = FALSE;
  if (abstau >= abscut && !ds->flag.planck && ds->nlyr > 1) {
    *lyrcut = TRUE;
  }
  if (!*lyrcut) {
    *ncut = ds->nlyr;
  }
  /*
   * Calculate Chapman function if spherical geometry, set expbea and ch for beam source.
   */
  if (ds->bc.fbeam > 0.) {
    CHTAU(0) = 0.;
    EXPBEA(0) = 1.;
    zenang    = acos(ds->bc.umu0)/DEG;

    if(ds->flag.spher == TRUE && ds->bc.umu0 < 0.) {
      EXPBEA(0) = exp(-c_chapman(1,0.,tauc,ds->nlyr,ds->zd,ds->dtauc,zenang,radius));
    }
    if (ds->flag.spher == TRUE) {
      for (lc = 1; lc <= *ncut; lc++) {
        taup        = TAUCPR(lc-1)+DTAUCPR(lc)/2.;
        CHTAU(lc  ) = c_chapman(lc, 0.0,      taucpr,ds->nlyr,ds->zd,dtaucpr,zenang,radius);
        chtau_tmp   = c_chapman(lc, 0.5,taucpr,ds->nlyr,ds->zd,dtaucpr,zenang,radius);
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
  /*
   * Set arrays defining location of user output levels within delta-m-scaled computational mesh
   */
  for (lu = 1; lu <= ds->ntau; lu++) {
    for (lc = 1; lc <= ds->nlyr-1; lc++) {
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
   * Set computational polar angle cosine for double gaussian
   * quadrature; cmu = 0.5, or  single gaussian quadrature; cmu = 1./sqrt(3
   * See KST for discussion of which is better for your specific applicatio
   */
  if(ds->flag.planck && ds->bc.fbeam == 0.) {
    *cmu = 0.5;
  }
  else {
    *cmu = sqrt(1./3.);
  }
  /*
   * Calculate planck functions
   */
  if (!ds->flag.planck) {
    *bplanck = 0.;
    *tplanck = 0.;
    /*
     * NOTE: If not using swappablecalloc() for dynamic memory allocation, need to zero-out
     *       pkag and pkagc here.
     */
  }
  else {
    *tplanck = emi_func(ds->wvnmlo,ds->wvnmhi,ds->bc.ttemp)*ds->bc.temis;
    *bplanck = emi_func(ds->wvnmlo,ds->wvnmhi,ds->bc.btemp);
    for (lev = 0; lev <= ds->nlyr; lev++) {
      PKAG(lev) = emi_func(ds->wvnmlo,ds->wvnmhi,TEMPER(lev));
    }
    for (lc = 1; lc <=ds->nlyr; lc++) {
      tempc     = .5*(TEMPER(lc-1)+TEMPER(lc));
      PKAGC(lc) = emi_func(ds->wvnmlo,ds->wvnmhi,tempc);
    }
  }

  return;
}

/*============================= end of c_twostr_set() ====================*/
