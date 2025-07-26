// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_disort_state_free() ===================*/

/*
 *  Free memory allocated by disort_state_alloc()
 */
void c_disort_state_free(disort_state *ds)
{
  if (ds->phi)    free(ds->phi);
  if (ds->umu)    free(ds->umu);
  if (ds->utau)   free(ds->utau);
  if (ds->temper) free(ds->temper);
  if (ds->pmom)   free(ds->pmom);
  if (ds->ssalb)  free(ds->ssalb);
  if (ds->dtauc)  free(ds->dtauc);
  if (ds->zd)     free(ds->zd);
  if (ds->flag.general_source == TRUE) {
    if (ds->gensrc)     free(ds->gensrc);
    if (ds->gensrcu)    free(ds->gensrcu);
  }
  if (!ds->flag.old_intensity_correction) {
    if (ds->nphase >= 1) {
      free(ds->mu_phase);
      free(ds->phase);
    }
  }

  switch(ds->flag.brdf_type) {
    case BRDF_RPV:
      free(ds->brdf.rpv);
    break;
#if HAVE_BRDF
    case BRDF_AMB:
      free(ds->brdf.ambrals);
    break;
    case BRDF_CAM:
      free(ds->brdf.cam);
    break;
#endif
    default:
      ;
    break;
  }

  return;
}

/*============================= end of c_disort_state_free() ============*/
