// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_disort_out_free() =====================*/

/*
 * Free memory allocated by disort_out_alloc()
 */
void c_disort_out_free(disort_state  *ds,
                       disort_output *out)
{

  if (out->trnmed) free(out->trnmed);
  if (out->albmed) free(out->albmed);
  if (out->u0u)    free(out->u0u);
  if (out->uu)     free(out->uu);
  if (out->rad)    free(out->rad);
  if ( ds->flag.output_uum )
    if (out->uum) free (out->uum);

  return;
}

/*============================= end of c_disort_out_free() ==============*/
