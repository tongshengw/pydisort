// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_inter() ===============================*/

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
 *
 *     Interpolates at the x-point arg from x-value array xarr and
 *     y-value array yarr. xarr and yarr are expected to have
 *     descending arguments, i.e. for atmospheric applications
 *     xarr typically holds the altitude and inter expects
 *     xarr(1) = altitude at top of atmosphere.
 *
 *     Input variables:
 *     dim       Array dimension of xarr and yarr
 *     npoints   No. points in arrays xarr and yarr
 *     itype     Interpolation type
 *     arg       Interpolation argument
 *     xarr      array of x values
 *     yarr      array of y values
 *
 *     Output variables:
 *     ynew      Interpolated function value at arg
 *     hh        gradient or scale height value
 *
 * This code was translated to c from fortran by Robert Buras
 *
 */

double c_inter( int     npoints,
		int     itype,
		double  arg,
		float  *xarr,
		double *yarr,
		double *hh )
{
  int iq=0, ip=0;
  double ynew=0.0;

  if ( arg <= XARR (1) && arg >= XARR (npoints) ) {
    for (iq=1;iq<=npoints-1;iq++)
      if ( arg <= XARR (iq) && arg >= XARR (iq+1) )
	ip=iq;
    if ( arg == XARR (npoints) )
      ip = npoints - 1;
  }
  else {
    if ( arg > XARR (1) )
      ip = 1;
    else {
      if ( arg < XARR (npoints) )
	ip = npoints - 1;
    }
  }

  /* Interpolate function value at arg from data points ip to ip+1 */

  switch(itype) {
  case 1:
    /*     exponential interpolation */
    if ( YARR (ip+1) == YARR (ip) ) {
      *hh = 0.0;
      ynew = YARR (ip);
    }
    else {
      *hh = -( XARR (ip+1) - XARR (ip) ) /
	log( YARR (ip+1) / YARR (ip));
      ynew = YARR (ip) * exp(- ( arg - XARR (ip) ) / *hh );
    }
    break;
  case 2:
    /*     linear interpolation */
    *hh = ( YARR (ip+1) - YARR (ip) ) / ( XARR (ip+1) - XARR (ip) );
    ynew = YARR (ip) + *hh * ( arg - XARR (ip) );
    break;
  default:
    printf("Error, unknown itype %d (line %d, function '%s' in '%s')\n",
	     itype, __LINE__, __func__, __FILE__);
    return -999.0;
  }

  return ynew;
}

/*============================= end of c_inter() ========================*/
