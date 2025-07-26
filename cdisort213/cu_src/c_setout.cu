// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_setout() ==============================*/

/*-------------------------------------------------------------------
 * Copyright (C) 1994 Arve Kylling
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 1, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * To obtain a copy of the GNU General Public License write to the
 * Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
 * USA.
 *-------------------------------------------------------------------
 * Linearly interpolate to get approximate tau corresponding to
 * altitude zout
 *
 * Input/output variables described in phodis.f
 *
 *
 * This code was translated to c from fortran by Robert Buras
 *
 */

int c_setout( float *sdtauc,
	      int    nlyr,
	      int    ntau,
	      float *sutau,
	      float *z,
	      float *zout )
{
  int itau=0, lc=0, itype=0;
  double hh=0.0;
  double *tauint=NULL;

  tauint = c_dbl_vector(0,nlyr+1,"tauint");

  if (tauint==NULL) {
    printf("Error allocating tauint!\n");
    return -1;
  }

  /* */

  TAUINT (1) = 0.0;
  for (lc=1; lc<=nlyr; lc++)
    TAUINT (lc+1) = TAUINT (lc) + SDTAUC (lc);

  itype = 2;

  for (itau=1; itau<=ntau; itau++)
    SUTAU (itau) = c_inter( nlyr+1, itype, ZOUT (itau),
			    z, tauint, &hh );

  free(tauint);

  return 0;
}

/*============================= end of c_setout() =======================*/
