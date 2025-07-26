// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_print_inputs() ==================*/

/*
 * Print values of twostream input variables
 */
void c_twostr_print_inputs(disort_state *ds,
                           int           deltam,
                           double       *flyr,
                           double       *gg,
                           int           lyrcut,
                           double       *oprim,
                           double       *tauc,
                           double       *taucpr)
{
  register int
    lu,lc;

  printf("\n\n"                 " ****************************************************************************************************\n"
                 " %s\n"
                 " ****************************************************************************************************\n",
                 ds->header);

  printf("\n No. streams = %4d     No. computational layers =%4d\n",ds->nstr,ds->nlyr);
  printf("%4d User optical depths :",ds->ntau);
  for (lu = 1; lu <= ds->ntau; lu++) {
    printf("%10.4f",UTAU(lu));
    if (lu%10 == 0) {
      printf("\n                          ");
    }
  }
  printf("\n");

  if (ds->flag.spher) {
    printf(" Pseudo-spherical geometry invoked\n");
  }

  if(!ds->flag.planck) {
    printf(" No thermal emission\n");
  }

  printf("    Incident beam with intensity =%11.3e and polar angle cosine = %8.5f\n"                 "    plus isotropic incident intensity =%11.3e\n",
                 ds->bc.fbeam,ds->bc.umu0,ds->bc.fisot);

  printf("    Bottom albedo (lambertian) =%8.4f\n",ds->bc.albedo);

  if(ds->flag.planck) {
    printf("    Thermal emission in wavenumber interval :%14.4f%14.4f\n"                   "    bottom temperature =%10.2f     top temperature =%10.2f    top emissivity =%8.4f\n",
                   ds->wvnmlo,ds->wvnmhi,ds->bc.btemp,ds->bc.ttemp,ds->bc.temis);
  }

  if(deltam) {
    printf(" Uses delta-m method\n");
  }
  else {
    printf(" Does not use delta-m method\n");
  }

  if(lyrcut) {
    printf(" Sets radiation = 0 below absorption optical depth 10\n");
  }

  if(ds->flag.planck) {
    printf("\n                                     <------------- delta-m --------------->"                   "\n                   total    single                           total    single"
                   "\n       optical   optical   scatter   truncated   optical   optical   scatter    asymm"
                   "\n         depth     depth    albedo    fraction     depth     depth    albedo   factor   temperature\n");
  }
  else {
    printf("\n                                     <------------- delta-m --------------->"                   "\n                   total    single                           total    single"
                   "\n       optical   optical   scatter   truncated   optical   optical   scatter    asymm"
                   "\n         depth     depth    albedo    fraction     depth     depth    albedo   factor\n");
  }

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if (ds->flag.planck) {
      printf("%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f%14.3f\n",                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),TAUCPR(lc)-TAUCPR(lc-1),TAUCPR(lc),OPRIM(lc),GG(lc),TEMPER(lc-1));
    }
    else {
      printf("%4d%10.4f%10.4f%10.5f%12.5f%10.4f%10.4f%10.5f%9.4f\n",                     lc,DTAUC(lc),TAUC(lc),SSALB(lc),FLYR(lc),TAUCPR(lc)-TAUCPR(lc-1),TAUCPR(lc),OPRIM(lc),GG(lc));
    }
  }

  if(ds->flag.planck) {
    printf("                                                                                     %14.3f\n",TEMPER(ds->nlyr));
  }

  return;
}

/*============================= end of c_twostr_print_inputs() ===========*/
