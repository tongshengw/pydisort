// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_xi_func() =============================*/

/*
   Calculates Xi function of eq. STWL (72)

         I N P U T   V A R I A B L E S

   umu1,2    cosine of zenith angle_1, _2
   tau       optical thickness of the layer

   NOTE: Original Fortran version also had argument umu3, but was only
         called for the case umu2 == umu3, so these two arguments are
         fused together here to reduce conditional testing.

   Called by- c_secondary_scat
 -------------------------------------------------------------------*/

double c_xi_func(double umu1,
               double umu2,
               double tau)
{
  double
    exp1,x1;

  x1   = (umu2-umu1)/(umu2*umu1);
  exp1 = exp(-tau/umu1);

  if (x1 != 0.) {
    return ((tau*x1-1.)*exp(-tau/umu2)+exp1)/(x1*x1*umu1*umu2);
  }
  else {
    return tau*tau*exp1/(2.*umu1*umu2);
  }
}

/*============================= end of c_xi_func() ======================*/
