// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_write_too_small_dim() ==================*/

/*
    Write name of too-small symbolic dimension and the value it should be
    increased to;  return TRUE

    Input :  quiet  = VERBOSE or QUIET
             dimnam = name of symbolic dimension which is too small
             minval = value to which that dimension should be increased
 ----------------------------------------------------------------------*/

int c_write_too_small_dim(int   quiet,
                          char const *dimnam,
                          int   minval)
{
  if (quiet != QUIET) {
    printf(" ****  Symbolic dimension %s should be increased to at least %d  ****\n",            dimnam,minval);
  }

  return TRUE;
}

/*============================= end of c_write_too_small_dim() ===========*/
