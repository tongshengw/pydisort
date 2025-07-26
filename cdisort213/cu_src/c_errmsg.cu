// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_errmsg() ===============================*/

/*
 * Print out a warning or error message;  abort if type == DS_ERROR
 */

#define MAX_WARNINGS 100

void c_errmsg(char const *messag,
              int   type)
{
  static int
    warning_limit = FALSE,
    num_warnings  = 0;

  if (type == DS_ERROR) {
    printf("\n ******* ERROR >>>>>>  %s\n",messag);
    __trap();
  }

  if (warning_limit) return;

  if (++num_warnings <= MAX_WARNINGS) {
    printf("\n ******* WARNING >>>>>>  %s\n",messag);
  }
  else {
    printf("\n\n >>>>>>  TOO MANY WARNING MESSAGES --  ','They will no longer be printed  <<<<<<<\n\n");
    warning_limit = TRUE;
  }

  return;
}

#undef MAX_WARNINGS

/*============================= end of c_errmsg() ========================*/
