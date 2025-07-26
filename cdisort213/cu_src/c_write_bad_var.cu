// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_write_bad_var() ========================*/

/*
   Write name of erroneous variable and return TRUE; count and abort
   if too many errors.

   Input : quiet  = VERBOSE or QUIET
           varnam = name of erroneous variable to be written
 ----------------------------------------------------------------------*/

int c_write_bad_var(int   quiet,
                    char const *varnam)
{
  const int
    maxmsg = 50;
  static int
    nummsg = 0;

  nummsg++;
  if (quiet != QUIET) {
    printf("\n ****  Input variable %s in error  ****\n",varnam);
    if (nummsg == maxmsg) {
      c_errmsg("Too many input errors.  Aborting...",DS_ERROR);
    }
  }

  return TRUE;
}

/*============================= end of c_write_bad_var() =================*/
