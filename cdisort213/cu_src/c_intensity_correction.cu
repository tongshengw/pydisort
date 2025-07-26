// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_intensity_correction() ================*/

/*
       Corrects intensity field by using Nakajima-Tanaka algorithm
       (1988). For more details, see Section 3.6 of STWL NASA report.
                I N P U T   V A R I A B L E S

       ds      Disort state variables
       dither  small multiple of machine precision
       flyr    separated fraction in delta-M method
       layru   index of UTAU in multi-layered system
       lyrcut  logical flag for truncation of computational layer
       ncut    total number of computational layers considered
       oprim   delta-M-scaled single-scatter albedo
       phirad  azimuthal angles in radians
       tauc    optical thickness at computational levels
       taucpr  delta-M-scaled optical thickness
       utaupr  delta-M-scaled version of UTAU

                O U T P U T   V A R I A B L E S

       out->UU  corrected intensity field; UU(IU,LU,J)
                 iu=1,ds->numu; lu=1,ds->ntau; j=1,ds->nphi

                I N T E R N A L   V A R I A B L E S

       ctheta  cosine of scattering angle
       dtheta  angle (degrees) to define aureole region as
                    direction of beam source +/- DTHETA
       phasa   actual (exact) phase function
       phasm   delta-M-scaled phase function
       phast   phase function used in TMS correction; actual phase
                    function divided by (1-FLYR*SSALB)
       pl      ordinary Legendre polynomial of degree l, P-sub-l
       plm1    ordinary Legendre polynomial of degree l-1, P-sub-(l-1)
       plm2    ordinary Legendre polynomial of degree l-2, P-sub-(l-2)
       theta0  incident zenith angle (degrees)
       thetap  emergent angle (degrees)
       ussndm  single-scattered intensity computed by using exact
                   phase function and scaled optical depth
                   (first term in STWL(68a))
       ussp    single-scattered intensity from delta-M method
                   (second term in STWL(68a))
       duims   intensity correction term from IMS method
                   (delta-I-sub-IMS in STWL(A.19))

   Called by- c_disort
   Calls- c_single_scat, c_secondary_scat
 -------------------------------------------------------------------*/

void c_intensity_correction(disort_state  *ds,
                            disort_output *out,
                            double         dither,
                            double        *flyr,
                            int           *layru,
                            int            lyrcut,
                            int            ncut,
                            double        *oprim,
                            double        *phasa,
                            double        *phast,
                            double        *phasm,
                            double        *phirad,
                            double        *tauc,
                            double        *taucpr,
                            double        *utaupr)
{
  register int
    iu,jp,k,lc,ltau,lu;
  double
    ctheta,dtheta,duims,pl,plm1,plm2,
    theta0=0,thetap=0,ussndm,ussp;

  dtheta = 10.;

  /*
   * Start loop over zenith angles
   */
  for (iu = 1; iu <= ds->numu; iu++) {
    if (UMU(iu) < 0.) {
      /*
       * Calculate zenith angles of incident and emerging directions
       */
      theta0 = acos(-ds->bc.umu0)/DEG;
      thetap = acos(UMU(iu))/DEG;
    }
    /*
     * Start loop over azimuth angles
     */
    for (jp = 1; jp <= ds->nphi; jp++) {
      /*
       * Calculate cosine of scattering angle, eq. STWL(4)
       */
      ctheta = -ds->bc.umu0*UMU(iu)+sqrt((1.-SQR(ds->bc.umu0))*(1.-SQR(UMU(iu))))*cos(PHIRAD(jp));
       /*
        * Initialize phase function
        */
      for (lc = 1; lc <= ncut; lc++) {
        PHASA(lc) = 1.;
        PHASM(lc) = 1.;
      }
      /*
       * Initialize Legendre poly. recurrence
       */

      plm1 = 1.;
      plm2 = 0.;
      for (k = 1; k <= ds->nmom; k++) {
        /*
         * Calculate Legendre polynomial of P-sub-l by upward recurrence
         */
        pl   = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
        plm2 = plm1;
        plm1 = pl;

        /*
         * Calculate actual phase function
         */
        for (lc = 1; lc <= ncut; lc++) {
          PHASA(lc) += (double)(2*k+1)*pl*PMOM(k,lc);
        }
        /*
         * Calculate delta-M transformed phase function
         */
        if (k <= ds->nstr-1) {
          for (lc = 1; lc <= ncut; lc++) {
            PHASM(lc) += (double)(2*k+1)*pl*(PMOM(k,lc)-FLYR(lc))/(1.-FLYR(lc));
          }
        }
      }
      /*
       * Apply TMS method, eq. STWL(68)
       */
      for (lc = 1; lc <= ncut; lc++) {
        PHAST(lc) = PHASA(lc)/(1.-FLYR(lc)*SSALB(lc));
      }
      for (lu = 1; lu <= ds->ntau; lu++) {
        if (!lyrcut || LAYRU(lu) < ncut) {
          ussndm        = c_single_scat(dither,LAYRU(lu),ncut,phast,ds->ssalb,taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          ussp          = c_single_scat(dither,LAYRU(lu),ncut,phasm,oprim,    taucpr,UMU(iu),ds->bc.umu0,UTAUPR(lu),ds->bc.fbeam);
          UU(iu,lu,jp) += ussndm-ussp;
        }
      }
      if (UMU(iu) < 0. && fabs(theta0-thetap) <= dtheta) {
        /*
         * Emerging direction is in the aureole (theta0 +/- dtheta).
         * Apply IMS method for correction of secondary scattering below top level.
         */
        ltau = 1;
        if (UTAU(1) <= dither) {
          ltau = 2;
        }
        for (lu = ltau; lu <= ds->ntau; lu++) {
          if(!lyrcut || LAYRU(lu) < ncut) {
            duims         = c_secondary_scat(ds,iu,lu,ctheta,flyr,LAYRU(lu),tauc);
	    UU(iu,lu,jp) -= duims;
          }
        }
      }
    } /* end loop over azimuth angles */
  } /* end loop over zenith angles */

  return;
}

/*============================= end of c_intensity_correction() =========*/
