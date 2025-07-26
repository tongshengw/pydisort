// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_interp_source() =======================*/

/*
    Interpolates source functions to user angles, eq. STWL(30)

    I N P U T      V A R I A B L E S:

       ds     :  Disort state variables
       cwt    :  Weights for Gauss quadrature over angle cosine
       delm0  :  Kronecker delta, delta-sub-m0
       gl     :  Delta-M scaled Legendre coefficients of phase function
                 (including factors 2l+1 and single-scatter albedo)
       mazim  :  Order of azimuthal component
       oprim  :  Single scattering albedo
       xr     :  Expansion of thermal source function, eq. STWL(24d); xr[].zero, xr[].one (see cdisort.h)
       ylm0   :  Normalized associated Legendre polynomial at the beam angle
       ylmc   :  Normalized associated Legendre polynomial at the quadrature angles
       ylmu   :  Normalized associated Legendre polynomial at the user angles
       zbs0   :  Solution vectors z-sub-zero of Eq. KS(10-11), used if pseudo-spherical
       zbs1   :  Solution vectors z-sub-one  of Eq. KS(10-11), used if pseudo-spherical
       zbsa   :  Alfa coefficient in Eq. KS(7), used if pseudo-spherical
       zee    :  Solution vectors Z-sub-zero, Z-sub-one of eq. SS(16), STWL(26a,b)
       zj     :  Solution vector Z-sub-zero after solving eq. SS(19), STWL(24b)
       zjg    :  Right-hand side vector  X-sub-zero in eq. KS(10), also the solution vector
                 Z-sub-zero after solving that system for a general source constant over a layer

    O U T P U T     V A R I A B L E S:

       zbeam  :  Incident-beam source function at user angles
       zu     :  Components 0 and 1 of a linear-in-optical-depth-dependent source (approximating the Planck emission source)
       zgu    :  General source function at user angles

   I N T E R N A L    V A R I A B L E S:

       psi  :   psi[].zero: Sum just after square bracket in eq. SD(9)
                psi[].one:  Sum in eq. STWL(31d)

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_interp_source(disort_state   *ds,
                     int             lc,
                     double         *cwt,
                     double          delm0,
                     double         *gl,
                     int             mazim,
                     double         *oprim,
                     double         *ylm0,
                     double         *ylmc,
                     double         *ylmu,
                     disort_pair    *psi,
                     disort_pair    *xr,
                     disort_pair    *zee,
                     double         *zj,
		     double         *zjg,
                     double         *zbeam,
		     disort_triplet *zbu,
		     disort_pair    *zbs,
		     double          zbsa,
		     double         *zgu,
                     disort_pair    *zu)
{
  register int
    iq,iu,jq;
  double
    fact,psum,psum0,psum1,sum,sum0,sum1;

  if (ds->bc.fbeam > 0.) {
    /*
     * Beam source terms; eq. SD(9)
     */
    if ( ds->flag.spher == TRUE ) {
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
	psum0 = 0.;
	psum1 = 0.;
	for (jq = 1; jq <= ds->nstr; jq++) {
	  psum0 +=  CWT(jq)*YLMC(iq,jq)*ZBS0(jq);
	  psum1 +=  CWT(jq)*YLMC(iq,jq)*ZBS1(jq);
	}
	PSI0(iq+1) = 0.5*GL(iq,lc)*psum0;
	PSI1(iq+1) = 0.5*GL(iq,lc)*psum1;
      }
      for (iu = 1; iu <= ds->numu; iu++) {
	sum0 = 0.;
	sum1 = 0.;
	for (iq = mazim; iq <= ds->nstr-1; iq++) {
	  sum0 += YLMU(iq,iu)*PSI0(iq+1);
	  sum1 += YLMU(iq,iu)*PSI1(iq+1);
	}
	ZB0U(iu,lc) = sum0 + ZB0U(iu,lc);
	ZB1U(iu,lc) = sum1 + ZB1U(iu,lc);
	ZBAU(iu,lc) = zbsa;
      }
    }
    else {
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
	psum = 0.;
	for (jq = 1; jq <= ds->nstr; jq++) {
	  psum += CWT(jq)*YLMC(iq,jq)*ZJ(jq);
	}
	PSI0(iq+1) = .5*GL(iq,lc)*psum;
      }
      fact = (2.-delm0)*ds->bc.fbeam/(4.*M_PI);
      for (iu = 1; iu <= ds->numu; iu++) {
	sum = 0.;
	for (iq = mazim; iq <= ds->nstr-1; iq++) {
	  sum += YLMU(iq,iu)*(PSI0(iq+1)+fact*GL(iq,lc)*YLM0(iq));
	}
	ZBEAM(iu,lc) = sum;
      }
    }
  }
  if (ds->flag.general_source > 0.) {
    /*
     * General source; eq. SD(9), KS(13)
     */
    for (iq = mazim; iq <= ds->nstr-1; iq++) {
      psum0 = 0.;
      for (jq = 1; jq <= ds->nstr; jq++) {
	psum0 +=  CWT(jq)*YLMC(iq,jq)*ZJG(jq);
      }
      PSI0(iq+1) = 0.5*GL(iq,lc)*psum0;
    }
    for (iu = 1; iu <= ds->numu; iu++) {
      sum0 = 0.;
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
	sum0 += YLMU(iq,iu)*PSI0(iq+1);
      }
      ZGU(iu,lc) = sum0 + GENSRCU(mazim,lc,iu);
    }
  }

  if (ds->flag.planck && mazim == 0) {
    /*
     * Thermal source terms, STWJ(27c), STWL(31c)
     */
    for (iq = mazim; iq <=ds->nstr-1; iq++) {
      psum0 = 0.;
      psum1 = 0.;
      for (jq = 1; jq <= ds->nstr; jq++) {
        psum0 += CWT(jq)*YLMC(iq,jq)*Z0(jq);
        psum1 += CWT(jq)*YLMC(iq,jq)*Z1(jq);
      }
      PSI0(iq+1) = .5*GL(iq,lc)*psum0;
      PSI1(iq+1) = .5*GL(iq,lc)*psum1;
    }
    for (iu = 1; iu <= ds->numu; iu++) {
      sum0 = 0.;
      sum1 = 0.;
      for (iq = mazim; iq <= ds->nstr-1; iq++) {
        sum0 += YLMU(iq,iu)*PSI0(iq+1);
        sum1 += YLMU(iq,iu)*PSI1(iq+1);
      }
      Z0U(iu,lc) = sum0+(1.-OPRIM(lc))*XR0(lc);
      Z1U(iu,lc) = sum1+(1.-OPRIM(lc))*XR1(lc);
    }
  }

  return;
}

/*============================= end of c_interp_source() ================*/
