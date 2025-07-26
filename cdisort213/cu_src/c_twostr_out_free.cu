// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_out_free() =====================*/

/*
 * Free memory allocated by twostr_out_alloc()
 */
void c_twostr_out_free(disort_state  *ds,
                       disort_output *out)
{
  if (out->rad) free(out->rad);

  return;
}

/*============================= end of c_twostr_out_free() ==============*/
