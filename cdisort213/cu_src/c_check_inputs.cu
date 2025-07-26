// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_check_inputs() ========================*/

/*
 * Checks the input dimensions and variables
 *
 * Calls- c_write_bad_var, c_dref, c_errmsg
 * Called by- c_disort
 */

int c_check_inputs(disort_state *ds,
		    int           scat_yes,
		    int           deltam,
		    int           corint,
		    double       *tauc,
		    int           callnum)
{
  int
    inperr = FALSE;
  register int
    irmu,iu,j,k,lc,lu, nu;
  double
    flxalb,rmu,umumin;

  if (ds->nstr < 2 || ds->nstr%2 != 0) {
    inperr = c_write_bad_var(VERBOSE,"ds.nstr");
  }
  if (ds->nstr == 2) {
    c_errmsg("check_inputs()--2 streams not recommended;\n\nUse specialized 2-stream code c_twostr() instead",DS_WARNING);
  }
  if (ds->nlyr < 1) {
    inperr = c_write_bad_var(VERBOSE,"ds.nlyr");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (DTAUC(lc) < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.dtauc");
    }
    if (SSALB(lc) < 0.0 || SSALB(lc) > 1.0) {
      inperr = c_write_bad_var(VERBOSE,"ds.ssalb");
    }
    if (ds->flag.ibcnd == GENERAL_BC) {
      if (ds->flag.planck) {
        if (lc == 1 && TEMPER(0) < 0.) {
          inperr = c_write_bad_var(VERBOSE,"ds.temper");
        }
        if (TEMPER(lc) < 0.) {
          inperr = c_write_bad_var(VERBOSE,"ds.temper");
        }
      }
    }
    else if (ds->flag.ibcnd == SPECIAL_BC) {
      ds->flag.planck = FALSE;
    }
    else {
      c_errmsg("check_inputs---unrecognized ds->flag.ibcnd",DS_ERROR);
    }
  }

  if (ds->nmom < 0 || (scat_yes  && ds->nmom < ds->nstr)) {
    inperr = c_write_bad_var(VERBOSE,"ds.nmom");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    for (k = 0; k <= ds->nmom; k++) {
      if (PMOM(k,lc) < -1. || PMOM(k,lc) > 1.) {
        inperr = c_write_bad_var(VERBOSE,"PMOM(k,lc)");
      }
    }
  }

  if( ds->flag.spher == TRUE ) {
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (ds->ZD(lc) > ds->ZD(lc-1)) {
        inperr     = c_write_bad_var(ds->flag.quiet,"zd");
      }
    }
  }

  if (ds->flag.ibcnd == GENERAL_BC) {
    if (ds->flag.usrtau) {
      if (ds->ntau < 1) {
        inperr = c_write_bad_var(VERBOSE,"ds.ntau");
      }
      for (lu = 1; lu <= ds->ntau; lu++) {
	/* Do a relative check to see if we are just beyond the bottom boundary */
	/* This might happen due to numerical rounding off problems.  ak20110224*/
        if (fabs(UTAU(lu)-TAUC(ds->nlyr)) <= 1.e-6*TAUC(ds->nlyr)) {
          UTAU(lu) = TAUC(ds->nlyr);
        }
        if(UTAU(lu) < 0. || UTAU(lu) > TAUC(ds->nlyr)) {
          inperr = c_write_bad_var(VERBOSE,"ds.utau");
        }
      }
    }
  }

  if (ds->flag.usrang) {
    if (ds->numu < 0) {
      inperr = c_write_bad_var(VERBOSE,"ds.numu");
    }
    if (!ds->flag.onlyfl && ds->numu == 0) {
      inperr = c_write_bad_var(VERBOSE,"ds.numu");
    }
    nu = ds->numu;
    if (ds->flag.ibcnd == SPECIAL_BC ) nu = ds->numu/2;
    for (iu = 1; iu <= nu; iu++) {
      if (UMU(iu) < -1. || UMU(iu) > 1. || UMU(iu) == 0.) {
        inperr = c_write_bad_var(VERBOSE,"ds.umu");
      }
      if (ds->flag.ibcnd == SPECIAL_BC && UMU(iu) < 0.) {
        inperr = c_write_bad_var(VERBOSE,"ds.umu");
      }
      if (iu > 1) {
        if (UMU(iu) < UMU(iu-1)) {
          inperr = c_write_bad_var(VERBOSE,"ds.umu");
        }
      }
    }
  }

  if (!ds->flag.onlyfl && ds->flag.ibcnd != SPECIAL_BC) {
    if (ds->nphi <= 0) {
      inperr = c_write_bad_var(VERBOSE,"ds.nphi");
    }
    for (j=1; j <=ds->nphi; j++) {
      if (PHI(j) < 0. || PHI(j) > 360.) {
        inperr = c_write_bad_var(VERBOSE,"ds.phi");
      }
    }
  }

  if (ds->flag.ibcnd != GENERAL_BC && ds->flag.ibcnd != SPECIAL_BC) {
    inperr = c_write_bad_var(VERBOSE,"ds.flag.ibcnd");
  }

  if (ds->flag.ibcnd == GENERAL_BC) {
    if (ds->bc.fbeam < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.fbeam");
    }
    else if (ds->bc.fbeam > 0.) {
      umumin = 0.;
      if( ds->flag.spher == TRUE ) {
	umumin = -1.;
      }
      if (ds->bc.umu0 <= umumin || ds->bc.umu0 > 1.) {
        inperr = c_write_bad_var(VERBOSE,"ds.bc.umu0");
      }
      if (ds->bc.phi0 < 0. || ds->bc.phi0 > 360.) {
        inperr = c_write_bad_var(VERBOSE,"ds.bc.phi0");
      }
    }

    if (ds->bc.fisot < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.fisot");
    }

    if (ds->flag.lamber) {
      if (ds->bc.albedo < 0. || ds->bc.albedo > 1.) {
        inperr = c_write_bad_var(VERBOSE,"ds.bc.albedo");
      }
    }
    else {
      /*
       * Make sure flux albedo at dense mesh of incident angles does not assume unphysical values
       */
      for (irmu = 0; irmu <= 100; irmu++) {
        rmu    = (double)irmu*0.01;
        flxalb = c_dref(ds->wvnmlo, ds->wvnmhi, rmu, ds->flag.brdf_type, &ds->brdf, callnum);
        if (flxalb < 0. || flxalb > 1.) {
          inperr = c_write_bad_var(VERBOSE,"bidir_reflectivity()");
        }
      }
    }
  }
  else {
    if (ds->bc.albedo < 0. || ds->bc.albedo > 1.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.albedo");
    }
  }

  if (ds->flag.planck && ds->flag.ibcnd != SPECIAL_BC) {
    if (ds->wvnmlo < 0. || ds->wvnmhi < ds->wvnmlo) {
      inperr = c_write_bad_var(VERBOSE,"ds.wvnmlo,hi");
    }
    if (ds->bc.temis < 0. || ds->bc.temis > 1.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.temis");
    }
    if (ds->bc.btemp < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.btemp");
    }
    if (ds->bc.ttemp < 0.) {
      inperr = c_write_bad_var(VERBOSE,"ds.bc.ttemp");
    }
  }

  if (ds->accur < 0. || ds->accur > 1.e-2) {
    inperr = c_write_bad_var(VERBOSE,"ds.accur");
  }

  if (inperr) {
    c_errmsg("DISORT--input and/or dimension errors",DS_WARNING);
    return 1;
  }

  if (ds->flag.planck && ds->flag.quiet == VERBOSE) {
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (fabs(TEMPER(lc)-TEMPER(lc-1)) > 10.) {
        c_errmsg("check_inputs--vertical temperature step may be too large for good accuracy",DS_WARNING);
      }
    }
  }
  if(!corint && (!ds->flag.onlyfl && ds->bc.fbeam > 0. && scat_yes && deltam)) {
    c_errmsg("check_inputs--intensity correction is off;\nintensities may be less accurate",DS_WARNING);
  }

  return 0;
}

/*============================= end of c_check_inputs() =================*/
