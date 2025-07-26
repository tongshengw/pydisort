// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_print_avg_intensities() ===============*/

/*
   Print azimuthally averaged intensities at user angles

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_print_avg_intensities(disort_state *ds,
			     disort_output *out)
{
  register int
    iu,iumax,iumin,
    lenfmt,lu,np,npass;

  if(ds->numu < 1) {
    return;
  }

  printf("\n\n *******  AZIMUTHALLY AVERAGED INTENSITIES (at user polar angles)  ********\n");
  lenfmt = 8;
  npass  = 1+(ds->numu-1)/lenfmt;
  printf("\n   Optical   Polar Angle Cosines"                 "\n     Depth");

  for (np = 1; np <= npass; np++) {
    iumin = 1+lenfmt*(np-1);
    iumax = IMIN(lenfmt*np,ds->numu);
    printf("\n          ");
    for (iu = iumin; iu <= iumax; iu++) {
      printf("%14.5f",UMU(iu));
    }
    printf("\n");

    for (lu = 1; lu <= ds->ntau; lu++) {
      printf("%10.4f",UTAU(lu));
      for (iu = iumin; iu <= iumax; iu++) {
        printf("%14.4e",U0U(iu,lu));
      }
      printf("\n");
    }
  }

  return;
}

/*============================= end of c_print_avg_intensities() ========*/
