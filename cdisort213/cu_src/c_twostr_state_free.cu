// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_state_free() ===================*/

/*
 *  Free memory allocated by twostr_state_alloc()
 */
void c_twostr_state_free(disort_state *ds)
{
  if (ds->utau)   free(ds->utau);
  if (ds->temper) free(ds->temper);
  if (ds->ssalb ) free(ds->ssalb);
  if (ds->dtauc ) free(ds->dtauc);
  if (ds->zd )    free(ds->zd);
  return;
}

/*============================= end of c_twostr_state_free() ============*/
