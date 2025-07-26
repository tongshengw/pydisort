// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_disort_state_alloc() ==================*/

/*
 * Dynamically allocate memory for disort input arrays, including
 * ones that the user can optionally ask disort() to calculate.
 */
void c_disort_state_alloc(disort_state *ds)
{
  int
    nu=0;

  ds->dtauc = c_dbl_vector(0,ds->nlyr,"ds->dtauc");
  ds->ssalb = c_dbl_vector(0,ds->nlyr,"ds->ssalb");
  /*
   * NOTE: PMOM is used in the code even when ds->nmom is not set by the user
   *       (such as when there is no scattering). Its first dimension needs to be
   *       at least 0:ds->nstr, hence we introduce ds->nmom_nstr in the C version.
   */
  ds->nmom_nstr = IMAX(ds->nmom,ds->nstr);
  ds->pmom      = c_dbl_vector(0,(ds->nmom_nstr+1)*ds->nlyr-1,"ds->pmom");

  if (ds->flag.ibcnd == SPECIAL_BC) {
    ds->flag.planck = FALSE;
    ds->flag.lamber = TRUE;
    ds->flag.usrtau = FALSE;
  }

  /* range 0 to nlyr */
  if (ds->flag.planck == TRUE) {
    ds->temper = c_dbl_vector(0,ds->nlyr,"ds->temper");
  }
  else {
    ds->temper = NULL;
  }

  if (ds->flag.general_source == TRUE) {
    ds->gensrc  = c_dbl_vector(0,ds->nstr*ds->nlyr*ds->nstr,"ds->gensrc");
    ds->gensrcu = c_dbl_vector(0,ds->nstr*ds->nlyr*ds->numu,"ds->gensrcu");
  }
  else {
    ds->gensrc  = NULL;
    ds->gensrcu = NULL;
  }

  if (ds->flag.usrtau == FALSE) {
    ds->ntau = ds->nlyr+1;
  }
  ds->utau = c_dbl_vector(0,ds->ntau-1,"ds->utau");

  //20130723ak Treat as in c_twostr_state_alloc. See comment there.
  //20130723ak Thanks to Tim for reporting this one.
  ds->zd        = c_dbl_vector(0,ds->nlyr+1,"ds->zd");

  /* range starts at 0 */
  nu = ds->numu;
  if ( (!ds->flag.usrang || ds->flag.onlyfl)) nu = ds->nstr;

  if (ds->flag.ibcnd == SPECIAL_BC)
    ds->umu = c_dbl_vector(0,2*nu,"ds->umu");
  else
    ds->umu = c_dbl_vector(0,nu,"ds->umu");

  if (ds->nphi >= 1) {
    ds->phi = c_dbl_vector(0,ds->nphi-1,"ds->nphi");
  }
  else {
    ds->phi = NULL;
  }

  if (!ds->flag.old_intensity_correction) {
    if (ds->nphase >= 1) {
      ds->mu_phase = c_dbl_vector(0,ds->nphase-1,"ds->mu_phase");
      ds->phase = c_dbl_vector(0,ds->nlyr*ds->nphase-1,"ds->phase");
    }
    else {
      ds->mu_phase = NULL;
      ds->phase = NULL;
    }
  }

  switch(ds->flag.brdf_type) {
    case BRDF_RPV:
      ds->brdf.rpv = (rpv_brdf_spec *)swappablecalloc(1,sizeof(rpv_brdf_spec));
      if (!ds->brdf.rpv) {
        c_errmsg("swappablecalloc error for ds->brdf.rpv",DS_ERROR);
      }
    break;
#if HAVE_BRDF
    case BRDF_AMB:
      ds->brdf.ambrals = (ambrals_brdf_spec *)swappablecalloc(1,sizeof(ambrals_brdf_spec));
      if (!ds->brdf.ambrals) {
        c_errmsg("swappablecalloc error for ds->brdf.ambrals",DS_ERROR);
      }
    break;
    case BRDF_CAM:
      ds->brdf.cam = (cam_brdf_spec *)swappablecalloc(1,sizeof(cam_brdf_spec));
      if (!ds->brdf.cam) {
        c_errmsg("swappablecalloc error for ds->brdf.cam",DS_ERROR);
      }
    break;
#endif
    default:
      ;
    break;
  }

  return;
}

/*============================= end of c_disort_state_alloc() ============*/
