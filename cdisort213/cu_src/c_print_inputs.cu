// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_print_inputs() ========================*/

/*
   Print values of input variables

   Called by- c_disort
 --------------------------------------------------------------------*/

void c_print_inputs(disort_state *ds,
                    double       *dtaucpr,
                    int           scat_yes,
                    int           deltam,
                    int           corint,
                    double       *flyr,
                    int           lyrcut,
                    double       *oprim,
                    double       *tauc,
                    double       *taucpr)
{
  register int
    iq,iu,j,k,lc,lu;

  printf("\n\n"                 " ****************************************************************************************************\n"
                 " DISORT: %s\n"
                 " ****************************************************************************************************\n",
                 ds->header);

  printf("\n No. streams =%4d     No. computational layers =%4d\n",ds->nstr,ds->nlyr);

  if (ds->flag.ibcnd != SPECIAL_BC) {
    printf("%4d User optical depths :",ds->ntau);
    for (lu = 1; lu <= ds->ntau; lu++) {
      printf("%10.4f",UTAU(lu));
      if (lu%10 == 0) {
        printf("\n                          ");
      }
    }
    printf("\n");
  }

  if (!ds->flag.onlyfl) {
    printf("%4d User polar angle cosines :",ds->numu);
    for (iu = 1; iu <= ds->numu; iu++) {
      printf("%9.5f",UMU(iu));
      if (iu%10 == 0) {
        printf("\n                               ");
      }
    }
    printf("\n");
  }

  if (!ds->flag.onlyfl && ds->flag.ibcnd != SPECIAL_BC) {
    printf("%4d User azimuthal angles :",ds->nphi);
    for (j = 1; j <= ds->nphi; j++) {
      printf("%9.2f",PHI(j));
      if (j%10 == 0) {
        printf("n                            ");
      }
    }
    printf("\n");
  }

  if (!ds->flag.planck || ds->flag.ibcnd == SPECIAL_BC) {
    printf(" No thermal emission\n");
  }

  if (ds->flag.spher == TRUE) {
    printf(" Pseudo-spherical geometry invoked\n");
  }

  if (ds->flag.general_source == TRUE) {
    printf(" Calculation with general source term\n");
  }

  if (ds->flag.ibcnd == GENERAL_BC) {
    printf(" Boundary condition flag: ds.flag.ibcnd = GENERAL_BC\n");
    printf("    Incident beam with intensity =%11.3e and polar angle cosine = %8.5f  and azimuth angle =%7.2f\n",                   ds->bc.fbeam,ds->bc.umu0,ds->bc.phi0);
    printf("    plus isotropic incident intensity =%11.3e\n",ds->bc.fisot);

    if (ds->bc.fluor > 0.0 ) {
      printf("    Bottom isotropic exiting intensity =%11.3e\n",ds->bc.fluor);
    }
    if (ds->flag.lamber) {
      printf("    Bottom albedo (Lambertian) =%8.4f\n",ds->bc.albedo);
    }
    else {
      printf("    Bidirectional reflectivity at bottom\n");
    }

    if(ds->flag.planck) {
      printf("    Thermal emission in wavenumber interval :%14.4f%14.4f\n",ds->wvnmlo,ds->wvnmhi);
      printf("    Bottom temperature =%10.2f    Top temperature =%10.2f    Top emissivity =%8.4f\n",                     ds->bc.btemp,ds->bc.ttemp,ds->bc.temis);
    }
  }
  else if (ds->flag.ibcnd == SPECIAL_BC) {
    printf(" Boundary condition flag: ds.flag.ibcnd = SPECIAL_BC\n");
    printf("    Isotropic illumination from top and bottom\n");
    printf("    Bottom albedo (Lambertian) =%8.4f\n",ds->bc.albedo);
  }
  else {
    c_errmsg("Unrecognized ds.flag.ibcnd",DS_WARNING);
  }

  if (deltam) {
    printf(" Uses delta-M method\n");
  }
  else {
    printf(" Does not use delta-M method\n");
  }

  if (corint) {
    printf(" Uses TMS/IMS method\n");
  }
  else {
    printf(" Does not use TMS/IMS method\n");
  }

  if (ds->flag.ibcnd == SPECIAL_BC) {
    printf(" Calculate albedo and transmissivity of medium vs. incident beam angle\n");
  }
  else if (ds->flag.onlyfl) {
    printf(" Calculate fluxes only\n");
  }
  else {
    printf(" Calculate fluxes and intensities\n");
  }

  printf(" Relative convergence criterion for azimuth series =%11.2e\n",ds->accur);

  if (lyrcut) {
    printf(" Sets radiation = 0 below absorption optical depth 10\n");
  }

  /*
   * Print layer variables (to read, skip every other line)
   */
  if(ds->flag.planck) {
    printf("\n                                     <------------- Delta-M --------------->");
    printf("\n                   Total    Single                           Total    Single");
    printf("\n       Optical   Optical   Scatter   Separated   Optical   Optical   Scatter    Asymm");
    printf("\n         Depth     Depth    Albedo    Fraction     Depth     Depth    Albedo   Factor   Temperature\n");
  }
  else {
    printf("\n                                     <------------- Delta-M --------------->");
    printf("\n                   Total    Single                           Total    Single");
    printf("\n       Optical   Optical   Scatter   Separated   Optical   Optical   Scatter    Asymm");
    printf("\n         Depth     Depth    Albedo    Fraction     Depth     Depth    Albedo   Factor\n");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {

    if (ds->flag.planck) {
      printf("%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f%14.3f\n",                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),DTAUCPR(lc),TAUCPR(lc),OPRIM(lc),PMOM(1,lc),TEMPER(lc-1));
    }
    else {
      printf("%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f\n",                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),DTAUCPR(lc),TAUCPR(lc),OPRIM(lc),PMOM(1,lc));
    }
  }
  if (ds->flag.planck) {
    printf("                                                                                     %14.3f\n",            TEMPER(ds->nlyr));
  }

  if (ds->flag.prnt[4] && scat_yes) {
    printf("\n Number of Phase Function Moments = %5d\n",ds->nmom+1);
    printf(" Layer   Phase Function Moments\n");
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (SSALB(lc) > 0.) {
        printf("%6d",lc);
        for (k = 0; k <= ds->nmom; k++) {
          printf("%11.6f",PMOM(k,lc));
          if ((k+1)%10 == 0) {
            printf("\n      ");
          }
        }
        printf("\n");
      }
    }
  }

  if (ds->flag.general_source == TRUE) {
    printf(" Calculation with general source term\n");
    j = 0;
    for (lc = 1; lc <= ds->nlyr; lc++) {
      printf("%4d%10.4f",lc,DTAUC(lc));
      for (iq = 1; iq <= ds->nstr; iq++) {
	printf("%13.6e",GENSRC(j,lc,iq));
      }
      printf("\n");
    }
  }


  return;
}

/*============================= end of c_print_inputs() =================*/
