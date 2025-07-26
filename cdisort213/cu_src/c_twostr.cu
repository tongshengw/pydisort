// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr() ===============================*/

/*
 Copyright (C) 1993, 1994, 1995 Arve Kylling

 C rewrite by Timothy E. Dowling (Univ. of Louisville)

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 1, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 To obtain a copy of the GNU General Public License write to the
 Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
 USA.

+---------------------------------------------------------------------+

     AUTHOR :  Arve Kylling (July 1993)
               Arve.Kylling@itek.norut.no

     REFERENCES (cited in the programs using the acronyms shown):

     DS: Dahlback, A. and K. Stamnes 1991: A new spherical
      model for computing the radiation field available
      for photolysis and heating at twilight, Planet.
      Space Sci. 39, 671-683.

     KS: Kylling, A., and K. Stamnes 1992: Efficient yet accurate
      solution of the linear transport  equation in the
      presence of internal sources: the exponential-linear
      approximation, J. Comp. Phys. 102, 265-276.

    KST: Kylling, A., K. Stamnes and S.-C. Tsay 1995: A reliable
      and efficient two-stream algorithm for radiative
      transfer; Documentation of accuracy in realistic
      layered media, in print, Journal of Atmospheric
      Chemistry 21, 115-150.

    STWJ: Stamnes, K., S.-C. Tsay, W. Wiscombe and K. Jayaweera
      1988: Numerically stable algorithm for discrete-
      ordinate-method radiative transfer in multiple
      scattering and emitting layered media, Appl.
      Optics., 27, 2502.

    WW: Wiscombe, W., 1977:  The Delta-M Method: Rapid Yet
      Accurate Radiative Flux Calculations, J. Atmos. Sci.
      34, 1408-1422

+---------------------------------------------------------------------+

    I n t r o d u c t o r y    n o t e

    (References are given as author-last-name strings, e.g., KST.)

    twostr() solves the radiative transfer equation in an absorbing,
    emitting and multiple scattering, layered pseudo-spherical
    medium in the two-stream approximation. For a discussion of the
    theory behind the present implementation see (KST).

    twostr() is based on the general n-stream algorithm DISORT
    described in Stamnes et al. (1988, STWJ), and incorporates
    all the advanced features of that algorithm. Furthermore it
    has been extended to include spherical geometry using the
    perturbation approach of Dahlback and Stamnes (1991). Relative
    to DISORT, it is both simplified and extended as follows:

     1) The only quantities calculated are mean intensities and fluxes.

     2) The medium may be taken to be pseudo-spherical (flag.spher is TRUE)

     3) Only Lambertian reflection at the bottom boundary is allowed


    General remarks about the structure of the input/output parameters

    The list of input variables is more easily comprehended if
    the following simple facts are borne in mind :

    * there is one vertical coordinate, measured in optical depth units;

    * the layers necessary for computational purposes are entirely
      decoupled from the levels at which the user desires results.

    The computational layering is usually constrained by the problem,
    in the sense that each computational layer must be reasonably
    homogeneous and not have a temperature variation of more than
    about 20 K across it (if thermal sources are considered).
    For example, a clear layer overlain by a cloud overlain by a
    dusty layer would suggest three computational layers.

    However the radiant quantities can be returned to the user at ANY
    level.  For example, the user may have picked 3 computational
    layers, but he can then request intensities from e.g. only the
    middle of the 2nd layer.

+---------------------------------------------------------------------+

    I n p u t    v a r i a b l e s

    Note on units:

       The radiant output units are determined by the sources of
    radiation driving the problem.  Lacking thermal emission, the
    radiant output units are the units of the sources ds.bc.fbeam and
    ds.bc.fisot.
       If thermal emission of any kind is included, subprogram planck_func2()
    determines the units.  The default planck_func2() has mks units [w/sq m].
    ds.bc.fbeam and ds.bc.fisot must have the same units as planck_func2() when
    thermal emission is present.


    ********  Computational layer structure  ********

        ===========================================================
        == Note:  Layers are numbered from the top boundary down ==
        ===========================================================

    ds.nlyr     Number of computational layers

    DTAUC(lc)   lc = 1 to ds.nlyr,
                optical depths of computational layers

    SSALB(lc)   lc = 1 to ds.nlyr,
                single-scatter albedos of computational layers

    GG(lc)      lc = 1 to ds.nlyr,
                asymmetry factor of computational layers
                Should be <= 1.0 (complete forward scattering) and
                >= -1.0 (complete backward scattering).
                NOTE. GG is changed by twostr() if deltam = TRUE.

    TEMPER(lev) lev = 0 to ds.nlyr, temperatures [K] of levels.
                (Note that temperature is specified at levels
                rather than for layers.)  Don't forget to put top
                temperature in 'TEMPER(0)', not 'TEMPER(1)'.  Top and
                bottom values do not need to agree with top and
                bottom boundary temperatures ds.bc.ttemp and ds.bc.btemp
                (i.e. slips are allowed).
                Needed only if ds.flag.planck is TRUE.

    ZD(lev)     lev = 0 to ds.nlyr, altitude of level above
                the ground, i.e. ZD(nlyr) = 0., the surface of
                the planet. Typically in units of (km)
                Must have same units as -radius-.
                Used to calculate the Chapman function when
                spherical geometry is needed.
                Needed only if flag.spher is TRUE.

    ds.wvnmlo,  Wavenumbers (inv cm) of spectral interval
      ds.wvnmhi ( used only for calculating Planck function )
                needed only if ds.flag.planck is true.
                If ds.wvnmlo < ds.wvnmhi the Planck function is
                integrated over this interval. If ds.wvnmlo ==  ds.wvnmhi
                the Planck function at wvnmlo is returned.


    ********  User level organization  ********

    ds.flag.usrtau = FALSE, radiant quantities are to be returned
                     at boundary of every computational layer.

                   = TRUE,  radiant quantities are to be returned
                     at user-specified optical depths, as follows:

    ds.ntau        Number of optical depths

    UTAU(lu)       lu = 1 to ds.ntau, user optical depths, in increasing order.
                   UTAU(ntau) must be no greater than the total optical depth of the medium.

     ******** Top and bottom boundary conditions  *********

    ds.bc.fbeam : Intensity of incident parallel beam at top boundary.
                  (units w/sq m if thermal sources active, otherwise
                  arbitrary units).  Corresponding incident flux
                  is  'umu0'  times 'fbeam'.  Note that this is an
                  infinitely wide beam, not a searchlight beam.

    ds.bc.umu0  : Polar angle cosine of incident beam.

    ds.bc.fisot : Intensity of top-boundary isotropic illumination.
                  (units w/sq m if thermal sources active, otherwise
                  arbitrary units).  Corresponding incident flux
                  is  pi (M_PI = 3.14159...)  times 'fisot'.

    ds.bc.albedo: Bottom-boundary albedo, bottom boundary is
                  assumed to be Lambert reflecting.

    ds.bc.btemp : Temperature of bottom boundary (K)  (bottom
                  emissivity is calculated from -albedo-,
                  so it need not be specified).
                  Needed only if -planck- is true.

    ds.bc.ttemp : Temperature of top boundary (K)
                  Needed only if -planck- is true.

    ds.bc.temis : Emissivity of top boundary
                  Needed only if -planck- is true.

    radius      : Distance from center of planet to the planets
                  surface (km)

    **********  Control flags  **************

    ds.flag.planck  = TRUE, include thermal emission
                      FALSE, ignore all thermal emission (saves computer time)
                     ( If ds.flag.planck = FALSE, it is not necessary to set any of
                      the variables having to do with thermal emission )

    ds.flag.prnt[0] = TRUE, print input variables
    ds.flag.prnt[1] = TRUE, print fluxes, mean intensities and flux divergence.

    deltam  = TRUE,  use delta-m method ( see Wiscombe, 1977 )
            = FALSE, don't use delta-m method
            In general intensities and fluxes will be more accurate
            for phase functions with a large forward peak (i.e.
            an asymmetry factor close to 1.) if 'deltam' is set true.

    ds.flag.spher = TRUE, spherical geometry accounted for. In this case
                    -radius- and -zd- must be set also. NOTE: this option
                    increases the execution time, hence use it only when
                    necessary if speed is of concern.
                  = FALSE, plane-parallel atmosphere assumed

    ds.header   : A 127- (or less) character header for prints


+---------------------------------------------------------------------+
               O u t p u t    v a r i a b l e s
+---------------------------------------------------------------------+

    == Note on units == If thermal sources are specified, fluxes come
                        out in [w/sq m] and intensities in [w/sq m/steradian].
                        Otherwise, the flux and intensity units are determined
                        by the units of -fbeam- and -fisot-.

    If ds.flag.usrtau = FALSE :

         ds.ntau      Number of optical depths at which radiant
                      quantities are evaluated ( = nlyr+1 )

         UTAU(lu)     lu = 1 to ntau, optical depths, in increasing
                      order, corresponding to boundaries of
                      computational layers (see -dtauc-)

    RFLDIR(lu)    :   Direct-beam flux (without delta-m scaling)

    RFLDN(lu)     :   Diffuse down-flux (total minus direct-beam)
                      (without delta-m scaling)

    FLUP(lu)      :   Diffuse up-flux

    DFDT(lu)      :   Flux divergence  d(net flux)/d(optical depth),
                      where 'net flux' includes the direct beam
                      (an exact result;  not from differencing fluxes)

    UAVG(lu)      :   Mean intensity (including the direct beam)

    IERROR(i)     :   Error flag array, if IERROR(i) is zero everything
                      is ok, otherwise twostr() found a fatal error, in this
                      case, twostr return immediately and reports the error in IERROR.

                      i =  1 : ds.nlyr <  1

                      i =  3 : dtauc   <  0.
                      i =  4 : ssalb   <  0. || ssalb > 1.
                      i =  5 : temper  <  0.
                      i =  6 : gg      < -1. || gg > 1.
                      i =  7 : ZD(lc)  >  ZD(lc-1)
                      i =  8 : ds.ntau <  1

                      i = 10 : UTAU(lu) < 0. || UTAU(lu) > TAUC(nlyr)

                      i = 12 : fbeam    < 0.
                      i = 13 : if flag.spher = FALSE
                                    umu0  < 0. || umu0 > 1.
                               if flag.spher = TRUE
                                    umu0  < 0. || umu0 > 1.
                      i = 14 : fisot   < 0.
                      i = 15 : albedo  < 0. || albedo > 1.
                      i = 16 : wvnmlo  < 0. || wvnmhi < wvnmlo
                      i = 17 : temis   < 0. || temis  > 1.
                      i = 18 : btemp   < 0.
                      i = 19 : ttemp   < 0.

                      i = 22 : !ds->flag.usrtau && ds->ntau < ds->nlyr+1

                     NOTE: i = 2, 9, 11, 20, and 21 are eliminated in the C version by the
                           change from static to dynamic memory allocation

+---------------------------------------------------------------------+

                 I/O variable specifications

+---------------------------------------------------------------------+
      Routines called (in order): c_twostr_check_inputs, c_twostr_set, c_twostr_print_inputs,
                                  c_twostr_solns, c_set_matrix, c_twostr_solve_bc, c_twostr_fluxes
+---------------------------------------------------------------------+

  Index conventions (for all loops and all variable descriptions):

     iq     :  For quadrature angles
     lu     :  For user levels
     lc     :  For computational layers (each having a different single-scatter albedo and/or phase function)
     lev    :  For computational levels
     ls     :  Runs from 0 to 2*ds->nlyr+1, ls = 1,2,3 refers to top, center and bottom of layer 1,
               ls = 3,4,5 refers to top, center and bottom of layer 2, etc.

+---------------------------------------------------------------------+

               I n t e r n a l    v a r i a b l e s

   B()...........Right-hand side vector of eqs. KST(38-41), set in twostr_solve_bc()
   bplanck.......Intensity emitted from bottom boundary
   CBAND().......Matrix of left-hand side of the linear system eqs. KST(38-41); in tridiagonal form
   CH(lc)........The Chapman-factor to correct for pseudo-spherical geometry in the direct beam.
   CHTAU(lc).....The optical depth in spherical geometry.
   cmu...........Computational polar angle, single or double Gaussian quadrature rule used, see twostr_set()
   EXPBEA(lc)....Transmission of direct beam in delta-m optical depth coordinates
   FLDIR(lu).....Direct beam flux (delta-m scaled); fl[].zero (see cdisort.h)
   FLDN(lu)......Diffuse down flux (delta-m scaled); fl[].one (see cdisort.h)
   FLYR(lc)......Truncated fraction in delta-m method
   KK(lc)........Eigenvalues in eq. KST(20)
   LAYRU(lu).....Computational layer in which user output level UTAU(lu) is located
   LL(iq,lc).....Constants of integration C-tilde in eqs. KST(42-43) obtained by solving eqs. KST(38-41)
   lyrcut........True, radiation is assumed zero below layer -ncut- because of almost complete absorption
   ncut..........Computational layer number in which absorption optical depth first exceeds abscut
   OPRIM(lc).....Single scattering albedo after delta-m scaling
   pass1.........TRUE on first entry, FALSE thereafter
   PKAG(0:lc)....Integrated Planck function for internal emission at layer boundaries
   PKAGC(lc).....Integrated Planck function for internal emission at layer center
   RR(lc)........Eigenvectors at polar quadrature angles.
   TAUC(0:lc)....Cumulative optical depth (un-delta-m-scaled)
   TAUCPR(0:lc)..Cumulative optical depth (delta-m-scaled if deltam = TRUE, otherwise equal to TAUC)
   tplanck.......Intensity emitted from top boundary
   U0C(iq,lu)....Azimuthally-averaged intensity
   UTAUPR(lu)....Optical depths of user output levels in delta-m coordinates;  equal to UTAU(lu) if no delta-m

   The following are members of the structure twostr_xyz:
   XB_0D(lc).....x-sub-zero-sup-minus in expansion of pseudo-spherical beam source, eq. KST(22)
   XB_0U(lc).....x-sub-zero-sup-plus  in expansion of pseudo-spherical beam source, eq. KST(22)
   XB_1D(lc).....x-sub-one-sup-minus  in expansion of pseudo-spherical beam source, eq. KST(22)
   XB_1U(lc).....x-sub-one-sup-plus   in expansion of pseudo-spherical beam source, eq. KST(22)
   XP_0(lc)......x-sub-zero in expansion of thermal source function; see eq. KST(22) (has no (mu) dependence)
   XP_1(lc)......x-sub-one  in expansion of thermal source function; see eq. KST(22) (has no (mu) dependence)
   YB_0D(lc).....y-sub-zero-sup-minus in eq. KST(23), solution for pseudo-spherical beam source
   YB_0U(lc).....y-sub-zero-sup-plus  in eq. KST(23), solution for pseudo-spherical beam source
   YB_1D(lc).....y-sub-one-sup-minus  in eq. KST(23), solution for pseudo-spherical beam source
   YB_1U(lc).....y-sub-one-sup-plus   in eq. KST(23), solution for pseudo-spherical beam source
   YP_0D(lc).....y-sub-zero-sup-minus in eq. KST(23), solution for thermal source
   YP_0U(lc).....y-sub-zero-sup-plus  in eq. KST(23), solution for thermal source
   YP_1D(lc).....y-sub-one-sup-minus  in eq. KST(23), solution for thermal source
   YP_1U(lc).....y-sub-one-sup-plus   in eq. KST(23), solution for thermal source
   ZB_A(lc)......Alpha coefficient in eq. KST(22) for pseudo-spherical beam source
   ZP_A(lc)......Alpha coefficient in eq. KST(22) for thermal source
*/

void c_twostr(disort_state  *ds,
              disort_output *out,
              int            deltam,
              double        *gg,
              int           *ierror,
              double         radius,
              emission_func_t emi_func)
{
  register int
    lc,ierr;
  int
    lyrcut,iret,ncut,nn;
  int
    *ipvt,
    *layru;
  ipvt  = (int *)swappablecalloc(ds->nstr*ds->nlyr, sizeof(int));
  layru = (int *)swappablecalloc(ds->ntau, sizeof(int));

  double
    cmu,bplanck,tplanck;
  double
    *b,*cband,*ch,*chtau,*dtaucpr,*expbea,*flyr,*ggprim,
    *kk,*ll,*oprim,*pkag,*pkagc,*rr,*tauc,*taucpr,*u0c,*utaupr;
  disort_pair
    *fl;
  twostr_xyz
    *ts;
  twostr_diag
    *diag;
  const double
    dither = 100.*DBL_EPSILON;

  /*
   * Allocate zeroed memory
   */
  b       = c_dbl_vector(0,ds->nstr*ds->nlyr-1,"b");
  cband   = c_dbl_vector(0,ds->nstr*ds->nlyr*(9*(ds->nstr/2)-2)-1,"cband");
  ch      = c_dbl_vector(0,ds->nlyr-1,"ch");
  chtau   = c_dbl_vector(0,(2*ds->nlyr+1)-1,"chtau");
  dtaucpr = c_dbl_vector(0,ds->nlyr-1,"dtaucpr");
  expbea  = c_dbl_vector(0,ds->nlyr,"expbea");
  flyr    = c_dbl_vector(0,ds->nlyr-1,"flyr");
  ggprim  = c_dbl_vector(0,ds->nlyr-1,"ggprim");
  kk      = c_dbl_vector(0,ds->nlyr-1,"kk");
  ll      = c_dbl_vector(0,ds->nlyr*ds->nstr-1,"ll");
  oprim   = c_dbl_vector(0,ds->nlyr-1,"oprim");
  pkag    = c_dbl_vector(0,ds->nlyr,"pkag");
  pkagc   = c_dbl_vector(0,ds->nlyr-1,"pkagc");
  rr      = c_dbl_vector(0,ds->nlyr-1,"rr");
  tauc    = c_dbl_vector(0,ds->nlyr,"tauc");
  taucpr  = c_dbl_vector(0,ds->nlyr,"taucpr");
  u0c     = c_dbl_vector(0,ds->ntau*ds->nstr-1,"u0c");
  utaupr  = c_dbl_vector(0,ds->ntau-1,"utaupr");
  /*
   * Using C structures to facilitate cache-aware memory allocation, which tends to
   * reduce cache misses and speed up computer execution.
   */
  fl   = (disort_pair *)swappablecalloc(ds->ntau,  sizeof(disort_pair)); if (!fl)   c_errmsg("twostr alloc error for fl",  DS_ERROR);
  ts   = (twostr_xyz  *)swappablecalloc(ds->nlyr,  sizeof(twostr_xyz )); if (!ts)   c_errmsg("twostr alloc error for ts",  DS_ERROR);
  diag = (twostr_diag *)swappablecalloc(2*ds->nlyr,sizeof(twostr_diag)); if (!diag) c_errmsg("twostr alloc error for diag",DS_ERROR);

  if(ds->flag.prnt[0]) {
    printf("\n\n\n\n"            " ************************************************************************************************************************\n"
            "                         Two stream method radiative transfer program, version 1.13\n"
            " ************************************************************************************************************************\n");
  }

  memset(ierror,0,TWOSTR_NERR*sizeof(int));

  /*
   * Calculate cumulative optical depth and dither single-scatter albedo to improve numerical behavior of
   * eigenvalue/vector computation
   */

  for (lc = 1; lc <= ds->nlyr; lc++) {
    if(SSALB(lc) == 1.) {
      SSALB(lc) = 1.-dither;
    }
    TAUC(lc) = TAUC(lc-1)+DTAUC(lc);
  }

  /*
   * Check input dimensions and variables
   */
  c_twostr_check_inputs(ds,gg,ierror,tauc);

  iret = 0;
  for (ierr = 1; ierr <= TWOSTR_NERR; ierr++) {
    if (IERROR(ierr) != 0) {
      iret = 1;
      if (ds->flag.quiet==VERBOSE) {
        printf("\ntwostr reports fatal error: %d\n",ierr);
      }
    }
  }
  if (iret == 1) {
    goto free_local_memory_and_return;
  }

 /*
  * Perform various setup operations
  */
  c_twostr_set(ds,&bplanck,ch,chtau,&cmu,deltam,dtaucpr,expbea,flyr,gg,ggprim,layru,&lyrcut,
             &ncut,&nn,oprim,pkag,pkagc,radius,tauc,taucpr,&tplanck,utaupr,emi_func);

  /*
   * Print input information
   */
  if (ds->flag.prnt[0]) {
    c_twostr_print_inputs(ds,deltam,flyr,gg,lyrcut,oprim,tauc,taucpr);
  }

  /*
   * Calculate the homogenous and particular solutions
   */
  c_twostr_solns(ds,ch,chtau,cmu,ncut,oprim,pkag,pkagc,taucpr,ggprim,kk,rr,ts);

  /*
   * Solve for constants of integration in homogeneous solution (general boundary conditions)
   */
  c_twostr_solve_bc(ds,ts,bplanck,cband,cmu,expbea,lyrcut,nn,ncut,tplanck,taucpr,kk,rr,ipvt,b,ll,diag);

  /*
   * Compute upward and downward fluxes, mean intensities and flux divergences.
   */
  c_twostr_fluxes(ds,ts,ch,cmu,kk,layru,ll,lyrcut,ncut,oprim,rr,taucpr,utaupr,out,u0c,fl);

  /*
   * Free allocated memory
   */
 free_local_memory_and_return:
  free(b),     free(cband), free(ch),  free(chtau),free(dtaucpr),free(expbea),free(flyr),free(fl);
  free(ggprim),free(kk),    free(ll),  free(oprim),free(pkag),   free(pkagc), free(rr),  free(tauc),free(taucpr);
  free(u0c),   free(utaupr),free(diag),free(ts);

  return;
}

/*============================= end of c_twostr() ========================*/
