// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_free_dbl_vector() =====================*/

/*
 * Frees memory allocated by dbl_vector().
 *
 * NOTE: If the array is zero-offset, can just use free().
 * NOTE: Argument nh is not used, but kept to match dbl_vector().
 */

void c_free_dbl_vector(double *m,
                       int     nl,
                       int     nh)
{
  int
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m      += nl_safe;
  free(m);

  return;
}

/*============================= end of c_free_dbl_vector() ==============*/
