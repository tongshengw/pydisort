// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_surface_bidir() =======================*/

/*
       Computes user's' surface bidirectional properties, STWL(41)

   I N P U T     V A R I A B L E S:

       ds     :  Disort input variables
       cmu    :  Computational polar angle cosines (Gaussian)
       delm0  :  Kronecker delta, delta-sub-m0
       mazim  :  Order of azimuthal component
       nn     :  Order of Double-Gauss quadrature (ds->nstr/2)
       callnum:  number of surface calls

    O U T P U T     V A R I A B L E S:

       bdr :  Fourier expansion coefficient of surface bidirectional
                 reflectivity (computational angles)
       rmu :  Surface bidirectional reflectivity (user angles)
       bem :  Surface directional emissivity (computational angles)
       emu :  Surface directional emissivity (user angles)

    I N T E R N A L     V A R I A B L E S:

       dref   :  Directional reflectivity
       gmu    :  The NMUG angle cosine quadrature points on (0,1)
                 NMUG is set in cdisort.h
       gwt    :  The NMUG angle cosine quadrature weights on (0,1)

   Called by- c_disort
   Calls- c_gaussian_quadrature, c_bidir_reflectivity
+---------------------------------------------------------------------*/

void c_surface_bidir(disort_state *ds,
                     double        delm0,
                     double       *cmu,
                     int           mazim,
                     int           nn,
                     double       *bdr,
                     double       *emu,
                     double       *bem,
                     double       *rmu,
		     int           callnum)
{
  static int
    pass1 = TRUE;
  register int
    iq,iu,jg,jq,k;
  double
    dref,sum;
  static double
    gmu[NMUG],gwt[NMUG];

  if (pass1) {
    pass1 = FALSE;
    c_gaussian_quadrature(NMUG/2,gmu,gwt);
    for (k = 1; k <= NMUG/2; k++) {
      GMU(k+NMUG/2) = -GMU(k);
      GWT(k+NMUG/2) =  GWT(k);
    }
  }

  memset(bdr,0,(ds->nstr/2)*((ds->nstr/2)+1)*sizeof(double));
  memset(bem,0,(ds->nstr/2)*sizeof(double));

  /*
   * Compute Fourier expansion coefficient of surface bidirectional reflectance
   * at computational angles eq. STWL (41)
   */
  if (ds->flag.lamber && mazim == 0) {
    for (iq = 1; iq <= nn; iq++) {
      BEM(iq) = 1.-ds->bc.albedo;
      for (jq = 0; jq <= nn; jq++) {
        BDR(iq,jq) = ds->bc.albedo;
      }
    }
  }
  else if (!ds->flag.lamber) {
    for (iq = 1; iq <= nn; iq++) {
      for (jq = 1; jq <= nn; jq++) {
        sum = 0.;
        for (k = 1; k <= NMUG; k++) {
          sum += GWT(k) *
	    c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, CMU(iq), CMU(jq),
				   M_PI * GMU(k), ds->flag.brdf_type, &ds->brdf, callnum)
	    * cos((double)mazim * M_PI * GMU(k) );
        }
        BDR(iq,jq) = .5*(2.-delm0)*sum;
      }
      if (ds->bc.fbeam > 0.) {
        sum = 0.;
        for(k = 1; k <= NMUG; k++) {
          sum += GWT (k) *
	    c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, CMU(iq), ds->bc.umu0,
				   M_PI * GMU(k), ds->flag.brdf_type, &ds->brdf, callnum )
	    * cos((double)mazim * M_PI * GMU(k) );
        }
        BDR(iq,0) = .5*(2.-delm0)*sum;
      }
    }
    if (mazim == 0) {
      /*
       * Integrate bidirectional reflectivity at reflection polar angle cosines -CMU- and incident angle
       * cosines -GMU- to get directional emissivity at computational angle cosines -CMU-.
       */
      for (iq = 1; iq <= nn; iq++) {
        dref = 0.;
        for (jg = 1; jg <= NMUG; jg++) {
          sum = 0.;
          for (k = 1; k <= NMUG/2; k++) {
            sum += GWT(k) * GMU(k) *
	      c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, CMU(iq), GMU(k),
				     M_PI * GMU(jg), ds->flag.brdf_type, &ds->brdf, callnum );
          }
          dref += GWT(jg)*sum;
        }
        BEM(iq) = 1.-dref;
      }
    }
  }
  /*
   * Compute Fourier expansion coefficient of surface bidirectional reflectance at user angles eq. STWL (41)
   */
  if(!ds->flag.onlyfl && ds->flag.usrang) {
    memset(emu,0,ds->numu*sizeof(double));
    memset(rmu,0,ds->numu*((ds->nstr/2)+1)*sizeof(double));
    for (iu = 1; iu <= ds->numu; iu++) {
      if (UMU(iu) > 0.) {
        if(ds->flag.lamber && mazim == 0) {
          for (iq = 0; iq <= nn; iq++) {
            RMU(iu,iq) = ds->bc.albedo;
          }
          EMU(iu) = 1.-ds->bc.albedo;
        }
        else if (!ds->flag.lamber) {
          for (iq = 1; iq <= nn; iq++) {
            sum = 0.;
            for (k = 1; k <= NMUG; k++) {
              sum += GWT(k) *
		c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, UMU(iu), CMU(iq),
				       M_PI * GMU(k), ds->flag.brdf_type, &ds->brdf, callnum )
		* cos( (double)mazim * M_PI * GMU(k) );
            }
            RMU(iu,iq) = .5*(2.-delm0)*sum;
          }
          if (ds->bc.fbeam > 0.) {
            sum = 0.;
            for (k = 1; k <= NMUG; k++) {
              sum += GWT(k) *
		c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, UMU(iu),
				       ds->bc.umu0, M_PI * GMU(k),
				       ds->flag.brdf_type, &ds->brdf, callnum )
		* cos( (double)mazim * M_PI * GMU(k) );
            }
            RMU(iu,0) = .5*(2.-delm0)*sum;
          }
          if (mazim == 0) {
            /*
             * Integrate bidirectional reflectivity at reflection angle cosines -UMU- and
             * incident angle cosines -GMU- to get directional emissivity at user angle cosines -UMU-.
             */
            dref = 0.;
            for (jg = 1; jg <= NMUG; jg++) {
              sum = 0.;
              for (k = 1; k <= NMUG/2; k++) {
                sum += GWT(k) * GMU(k) *
		  c_bidir_reflectivity ( ds->wvnmlo, ds->wvnmhi, UMU(iu), GMU(k),
					 M_PI*GMU(jg), ds->flag.brdf_type, &ds->brdf, callnum );
              }
              dref += GWT(jg)*sum;
            }
            EMU(iu) = 1.-dref;
          }
        }
      }
    }
  }

  return;
}

/*============================= end of c_surface_bidir() ================*/
