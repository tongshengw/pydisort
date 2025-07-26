// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_check_inputs() ==================*/

/*
 * Checks the twostr input dimensions and variables
 */

void c_twostr_check_inputs(disort_state *ds,
                           double       *gg,
                           int          *ierror,
                           double       *tauc)
{
  int
    inperr,lc,lu;
  double
    umumin;

  inperr = FALSE;

  if (ds->nlyr < 1) {
    inperr    = c_write_bad_var(ds->flag.quiet,"nlyr");
    IERROR(1) = 1;
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (DTAUC(lc) < 0.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"dtauc");
      IERROR(3) += 1;
    }
    if (SSALB(lc) < 0. || SSALB(lc) > 1.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"ssalb");
      IERROR(4) += 1;
    }
    if (ds->flag.planck) {
      if (lc == 1 && TEMPER(0) < 0.) {
        inperr     = c_write_bad_var(ds->flag.quiet,"temper");
        IERROR(5) += 1;
      }
      if (TEMPER(lc) < 0.) {
        inperr     = c_write_bad_var(ds->flag.quiet,"temper");
        IERROR(5) += 1;
      }
    }
    if (GG(lc) < -1. || GG(lc) > 1.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"gg");
      IERROR(6) += 1;
    }
  }

  if(ds->flag.spher==TRUE) {
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (ds->ZD(lc) > ds->ZD(lc-1)) {
        inperr     = c_write_bad_var(ds->flag.quiet,"zd");
        IERROR(7) += 1;
      }
    }
  }

  if (ds->flag.usrtau) {
    if (ds->ntau < 1) {
      inperr    = c_write_bad_var(ds->flag.quiet,"ntau");
      IERROR(8) = 1;
    }
    for (lu = 1; lu <= ds->ntau; lu++) {
      if (fabs(UTAU(lu)-TAUC(ds->nlyr)) <= 1.e-6*TAUC(ds->nlyr)) { /* relative check copied from c_check_inputs() */
        UTAU(lu)= TAUC(ds->nlyr);
      }
      if (UTAU(lu) < 0. || UTAU(lu) > TAUC(ds->nlyr)) {
        inperr      = c_write_bad_var(ds->flag.quiet,"utau");
        IERROR(10) += 1;
      }
    }
  }

  if (ds->bc.fbeam < 0.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"fbeam");
    IERROR(12) = 1;
  }

  umumin = 0.;
  if(ds->flag.spher==TRUE) {
    umumin = -1.;
  }

  if (ds->bc.fbeam > 0. && (ds->bc.umu0 <= umumin || ds->bc.umu0 > 1.)) {
    inperr     = c_write_bad_var(ds->flag.quiet,"umu0");
    IERROR(13) = 1;
  }
  if (ds->bc.fisot < 0.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"fisot");
    IERROR(14) = 1;
  }
  if (ds->bc.albedo < 0. || ds->bc.albedo > 1.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"albedo");
    IERROR(15) = 1;
  }

  if(ds->flag.planck) {
    if (ds->wvnmlo < 0. || ds->wvnmhi < ds->wvnmlo) {
      inperr     = c_write_bad_var(ds->flag.quiet,"wvnmlo,hi");
      IERROR(16) = 1;
    }
    if (ds->bc.temis < 0. || ds->bc.temis > 1.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"temis");
      IERROR(17) = 1;
    }
    if (ds->bc.btemp < 0.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"btemp");
      IERROR(18) = 1;
    }
    if (ds->bc.ttemp < 0.) {
      inperr     = c_write_bad_var(ds->flag.quiet,"ttemp");
      IERROR(19) = 1;
    }
  }

  if (!ds->flag.usrtau && ds->ntau < ds->nlyr+1) {
    inperr = c_write_too_small_dim(ds->flag.quiet,"ds.ntau",ds->nlyr+1);
    IERROR(22) = 1;
  }

  if (ds->bc.fluor < 0.) {
    inperr     = c_write_bad_var(ds->flag.quiet,"fluor");
    IERROR(23) = 1;
  }

  if (inperr) {
    c_errmsg("twostr_check_inputs--input and/or dimension errors",DS_ERROR);
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (ds->flag.planck && fabs(TEMPER(lc)-TEMPER(lc-1)) > 50. && ds->flag.quiet==VERBOSE) {
      c_errmsg("twostr_check_inputs--vertical temperature step may be too large for good accuracy",DS_WARNING);
    }
  }

  return;
}

/*============================= end of c_twostr_check_inputs() ===========*/
