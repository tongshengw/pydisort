// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_out_alloc() ====================*/

/*
 *   Dynamically allocate memory for twostr output arrays
 */
void c_twostr_out_alloc(disort_state  *ds,
                        disort_output *out)
{
  out->rad = (disort_radiant *)swappablecalloc(ds->ntau,sizeof(disort_radiant));
  if (!out->rad) {
    c_errmsg("disort_out_alloc---error allocating out->rad array",DS_ERROR);
  }

  return;
}

/*============================= end of c_twostr_out_alloc() =============*/
