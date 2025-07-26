// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_gaussian_quadrature_test() ============*/

int c_gaussian_quadrature_test(int nstr, float *sza, double umu0)
{

  /* Test if the solar zenith angle coincides with one of
     the computational angles */

  int nn=0, iq=0, result=0;
  double umu0s=0.0, *cmu=NULL, *cwt=NULL;

  cmu = c_dbl_vector(0,nstr,"cmu");
  if (cmu==NULL) {
    printf("Error allocating cmu!\n");
    return -1;
  }

  cwt = c_dbl_vector(0,nstr,"cwt");
  if (cwt==NULL) {
    printf("Error allocating cwt!\n");
    return -1;
  }

  nn = nstr / 2.0;

  c_gaussian_quadrature ( nn, cmu, cwt );

  for (iq=1; iq<=nn; iq++) {
    if( fabs( (umu0 - CMU (iq)) / umu0 ) < 1.0e-4 ) {
      umu0s = umu0;
      if ( umu0 < CMU (iq) )
	umu0  = CMU (iq) * (1. - 1.1e-4);
      else
	umu0  = CMU (iq) * (1. + 1.1e-4);

      *sza   = acos (umu0)/DEG;
      printf("%s %s %s %f %s %f\n",	      "******* WARNING >>>>>> \n",
	      "SETDIS--beam angle=computational angle;\n",
	      "******* changing cosine of solar zenith angle, umu0, from ",
	      umu0s, "to", umu0 );
      result=-1;
    }
  }

  free(cwt);
  free(cmu);
  return result;
}

/*============================= end of c_gaussian_quadrature_test() =====*/
