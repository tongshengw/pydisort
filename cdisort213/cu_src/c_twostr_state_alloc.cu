// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_state_alloc() ==================*/

/*
 * Dynamically allocate memory for twostr input arrays.
 */
void c_twostr_state_alloc(disort_state *ds)
{
  /* Set to two streams */
  ds->nstr = 2;

  /* Set flags not controlled by user */
  ds->flag.prnt[2] = FALSE;
  ds->flag.prnt[3] = FALSE;
  ds->flag.prnt[4] = FALSE;
  ds->flag.onlyfl  = TRUE;

  ds->dtauc = c_dbl_vector(0,ds->nlyr-1,"ds->dtauc");
  ds->ssalb = c_dbl_vector(0,ds->nlyr-1,"ds->ssalb");

  /* range 0 to nlyr */
  if (ds->flag.planck == TRUE) {
    ds->temper = c_dbl_vector(0,ds->nlyr,"ds->temper");
  }
  else {
    ds->temper = NULL;
  }

  if (ds->flag.usrtau == FALSE) {
    ds->ntau = ds->nlyr+1;
  }
  ds->utau = c_dbl_vector(0,ds->ntau-1,"ds->utau");

  //20120820ak Tim says: if spher is false
  //20120820ak during allocation, it has a seg fault later if you turn spher
  //20120820ak on and try to use that functionality.  Probably would be
  //20120820ak better to just allocate this little guy regardless of the status of spher.
  //20120820ak So I commented the following.
  //20120820ak  if (ds->flag.spher == TRUE) {
  ds->zd        = c_dbl_vector(0,ds->nlyr+1,"ds->zd");
  //20120820ak  }

  return;
}

/*============================= end of c_twostr_state_alloc() ============*/
