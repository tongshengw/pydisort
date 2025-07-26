// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
/*============================= c_print_albtrans() =======================*/

/*
   Print planar albedo and transmissivity of medium as a function of
   incident beam angle

   Called by- c_albtrans
 --------------------------------------------------------------------*/

void c_print_albtrans(disort_state  *ds,
                      disort_output *out)
{
  register int
    iu;

  printf("\n\n\n *******  Flux Albedo and/or Transmissivity of entire medium  ********\n");
  printf(" Beam Zen Ang   cos(Beam Zen Ang)      Albedo   Transmissivity\n");
  for (iu = 1; iu <= ds->numu; iu++) {
    printf("%13.4f%20.6f%12.5f%17.4e\n",acos(UMU(iu))/DEG,UMU(iu),ALBMED(iu),TRNMED(iu));
  }

  return;
}

/*============================= end of c_print_albtrans() ================*/
