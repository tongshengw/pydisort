// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_int_vector() ==========================*/

/*
 * Allocates memory for a 1D integer array with range [nl..nh].
 *
 * NOTE: swappablecalloc() zeros the memory it allocates.
 */

int *c_int_vector(int  nl,
		  int  nh,
		  char const *name)
{
  unsigned int
    len_safe;
  int
    nl_safe, nh_safe;
  int
    *m;

  if (nh < nl) {
    printf("\n\n**error:%s, variable %s, range (%d,%d)\n","int_vector",name,nl,nh);
    __trap();
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (int *)swappablecalloc(len_safe,sizeof(int));

  if (!m) {
    c_errmsg("int_vector---alloc error",DS_ERROR);
  }
  m -= nl_safe;

  return m;
}

/*============================= end of c_int_vector() ===================*/
