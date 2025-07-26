// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_print_intensities() ===================*/

/*
   Prints the intensity at user polar and azimuthal angles
   All arguments are disort state or output variables

   Called by- c_disort
 -------------------------------------------------------------------*/

void c_print_intensities(disort_state  *ds,
                         disort_output *out)
{
  register int
    iu,j,jmax,jmin,lenfmt,lu,np,npass;

  if (ds->nphi < 1) {
    return;
  }

  printf("\n\n *********  I N T E N S I T I E S  *********\n");
  lenfmt = 10;
  npass  = 1+(ds->nphi-1)/lenfmt;
  printf("\n             Polar   Azimuth angles (degrees)");
  printf("\n   Optical   Angle");
  printf("\n    Depth   Cosine\n");
  for (lu = 1; lu <= ds->ntau; lu++) {
    for (np = 1; np <= npass; np++) {
      jmin = 1+lenfmt*(np-1);
      jmax = IMIN(lenfmt*np,ds->nphi);
      printf("\n                  ");
      for (j = jmin; j <= jmax; j++) {
        printf("%11.2f",PHI(j));
      }
      printf("\n");
      if (np == 1) {
        printf("%10.4f%8.4f",UTAU(lu),UMU(1));
        for (j = jmin; j <= jmax; j++) {
          printf("%11.3e",UU(1,lu,j));
        }
        printf("\n");
      }
      else {
        printf("          %8.4f",UMU(1));
        for (j = jmin; j <= jmax; j++) {
          printf("%11.3e",UU(1,lu,j));
        }
        printf("\n");
      }
      for (iu = 2; iu <= ds->numu; iu++) {
        printf("          %8.4f",UMU(iu));
        for (j = jmin; j <= jmax; j++) {
          printf("%11.3e",UU(iu,lu,j));
        }
        printf("\n");
      }
    }
  }

  return;
}

/*============================= end of c_print_intensities() ============*/
