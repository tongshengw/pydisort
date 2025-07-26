// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_disort_out_alloc() ====================*/

/*
 *   Dynamically allocate memory for disort output arrays
 */
void c_disort_out_alloc(disort_state  *ds,
                        disort_output *out)
{

  int
    nu;

  out->rad = (disort_radiant *)swappablecalloc(ds->ntau,sizeof(disort_radiant));

  if (!out->rad) {
    c_errmsg("disort_out_alloc---error allocating out->rad array",DS_ERROR);
  }
  nu = ds->numu;
  if ( (!ds->flag.usrang || ds->flag.onlyfl)) {
    nu = ds->nstr;
  }
  out->uu = c_dbl_vector(0,ds->nphi*nu*ds->ntau,"out->uu");

  out->u0u = c_dbl_vector(0,ds->ntau*nu,"out->u0u");

  if ( ds->flag.output_uum )
    out->uum = c_dbl_vector(0,ds->nstr*nu*ds->ntau,"out->uum");

  if (ds->flag.ibcnd == SPECIAL_BC) {
    out->albmed = c_dbl_vector(0,ds->numu,"out->albmed");
    out->trnmed = c_dbl_vector(0,ds->numu,"out->trnmed");
  }
  else {
    out->albmed = NULL;
    out->trnmed = NULL;
  }

  return;
}

/*============================= end of c_disort_out_alloc() =============*/
