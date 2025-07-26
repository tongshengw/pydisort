// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_secondary_scat() ======================*/

/*
   Calculates secondary scattered intensity of eq. STWL (A7)

                I N P U T   V A R I A B L E S

        ds      Disort state variables
        iu      index of user polar angle
        lu      index of user level
        ctheta  cosine of scattering angle
        flyr    separated fraction f in Delta-M method
        layru   index of utau in multi-layered system
        tauc    cumulative optical depth at computational layers

                I N T E R N A L   V A R I A B L E S

        pspike  2*p"-p"*p", where p" is the residual phase function
        wbar    mean value of single scattering albedo
        fbar    mean value of separated fraction f
        dtau    layer optical depth
        stau    sum of layer optical depths between top of atmopshere and layer layru

   Called by- c_intensity_correction
   Calls- c_xi_func
 -------------------------------------------------------------------*/

double c_secondary_scat(disort_state *ds,
                        int           iu,
                        int           lu,
                        double        ctheta,
                        double       *flyr,
                        int           layru,
                        double       *tauc)
{
  register int
    k,lyr;
  const double
    tiny = 1.e-4;
  double
    dtau,fbar,gbar,pl,plm1,plm2,pspike,
    stau,umu0p,wbar;
  register double
    tmp;

  /*
   * Calculate vertically averaged value of single scattering albedo and separated
   * fraction f, eq. STWL (A.15)
   */
  dtau = UTAU(lu)-TAUC(layru-1);
  wbar = SSALB(layru)*dtau;
  fbar = FLYR(layru)*wbar;
  stau = dtau;
  for (lyr = 1; lyr <= layru-1; lyr++) {
    wbar += DTAUC(lyr)*SSALB(lyr);
    fbar += DTAUC(lyr)*SSALB(lyr)*FLYR(lyr);
    stau += DTAUC(lyr);
  }

  if (wbar <= tiny || fbar <= tiny || stau <= tiny || ds->bc.fbeam <= tiny) {
    return 0.;
  }

  fbar /= wbar;
  wbar /= stau;
  /*
   * Calculate pspike = (2p"-p"*p")
   */
  pspike = 1.;
  gbar   = 1.;
  plm1   = 1.;
  plm2   = 0.;
  /*
   * pspike for l <= 2n-1
   */
  for (k = 1; k <= ds->nstr-1; k++) {
    pl      = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
    plm2    = plm1;
    plm1    = pl;
    pspike += gbar*(2.-gbar)*(double)(2*k+1)*pl;
  }
  /*
   * pspike for l > 2n-1
   */
  for (k = ds->nstr; k <= ds->nmom; k++) {
    pl   = ((double)(2*k-1)*ctheta*plm1-(double)(k-1)*plm2)/k;
    plm2 = plm1;
    plm1 = pl;
    dtau = UTAU(lu)-TAUC(layru-1);
    gbar = PMOM(k,layru)*SSALB(layru)*dtau;
    for (lyr = 1; lyr <= layru-1; lyr++) {
      gbar += PMOM(k,lyr)*SSALB(lyr)*DTAUC(lyr);
    }
    tmp = fbar*wbar*stau;
    if (tmp <= tiny) {
      gbar = 0.;
    }
    else {
      gbar /= tmp;
    }
    pspike += gbar*(2.-gbar)*(double)(2*k+1)*pl;
  }
  umu0p = ds->bc.umu0/(1.-fbar*wbar);
  /*
   * Calculate IMS correction term, eq. STWL (A.13)
   */
  return ds->bc.fbeam/(4.*M_PI)*SQR(fbar*wbar)/(1.-fbar*wbar)*pspike*c_xi_func(-UMU(iu),umu0p,UTAU(lu));
}

/*============================= end of c_secondary_scat() ===============*/
