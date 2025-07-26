// includes
#include<alloc.h>
#include<cdisort.h>
#include<locate.h>

DISPATCH_MACRO
#undef  KK
#define KK(lyu) kk[lyu-1]
/*============================= c_twostr_solve_bc() ======================*/

/*
 Construct right-hand side vector -b- for general boundary conditions
 and solve system of equations obtained from the boundary conditions
 and the continuity-of-intensity-at-layer-interface equations.

 Routines called: c_sgbfa, c_sgbsl

 I n p u t      v a r i a b l e s:

       ds       : 'Disort' state variables
       ts       :  twostr_xyz structure variables (see cdisort.h)
       bplanck  :  Bottom boundary thermal emission
       cband    :  Left-hand side matrix of linear system eqs. KST(38-41)
                   in banded form required by linpack solution routines
       cmu      :  Abscissa for gauss quadrature over angle cosine
       expbea   :  Transmission of incident beam, EXP(-taucpr/ch)
       lyrcut   :  Logical flag for truncation of comput. layer
       ncol     :  Counts of columns in -cband-
       nn       :  Order of double-gauss quadrature (nstr/2)
       ncut     :  Total number of computational layers considered
       tplanck  :  Top boundary thermal emission
       taucpr   :  Cumulative optical depth (delta-m-scaled)
       kk       :
       rr       :
       ipvt     :

 O u t p u t     v a r i a b l e s:

       b        :  Right-hand side vector of eqs. KST(38-41) going into
                   sgbsl; returns as solution vector of eqs. KST(38-41)
                   constants of integration
       ll       :  Permanent storage for -b-, but re-ordered

 I n t e r n a l    v a r i a b l e s:

       diag     : diag[].super, diag[].on, diag[].sub

 ---------------------------------------------------------------------*/

void c_twostr_solve_bc(disort_state *ds,
                       twostr_xyz   *ts,
                       double        bplanck,
                       double       *cband,
                       double        cmu,
                       double       *expbea,
                       int           lyrcut,
                       int           nn,
                       int           ncut,
                       double        tplanck,
                       double       *taucpr,
                       double       *kk,
                       double       *rr,
                       int          *ipvt,
                       double       *b,
                       double       *ll,
                       twostr_diag  *diag)
{
  int
    info;
  register int
    irow,lc,nloop,nrow,job;
  double
    wk0,wk1,wk,rpp1_m,rp_m,rpp1_p,rp_p,sum,refflx;
  register double
    fact1,fact2,fact3,fact4;

  /*
   * First top row, top boundary condition
   */
  irow = 1;
  lc   = 1;
  /*
   * SUBD(irow) is undefined
   */
  DIAG(irow)   = RR(lc)*exp(-KK(lc)*TAUCPR(lc));
  SUPERD(irow) = 1.;
  /*
   * next from layer no. 2 to nlyr-1
   */
  nloop = ncut-1;
  for (lc = 1; lc <= nloop; lc++) {
    irow++;
    wk0          = exp(-KK(lc  )*(TAUCPR(lc  )-TAUCPR(lc-1)));
    wk1          = exp(-KK(lc+1)*(TAUCPR(lc+1)-TAUCPR(lc  )));
    SUBD(irow)   = 1.-RR(lc)*RR(lc+1);
    DIAG(irow)   = (RR(lc)-RR(lc+1))*wk0;
    SUPERD(irow) = -(1.-SQR(RR(lc+1)))*wk1;
    irow++;
    SUBD(irow)   = (1.-SQR(RR(lc)))*wk0;
    DIAG(irow)   = (RR(lc)-RR(lc+1))*wk1;
    SUPERD(irow) = -(1.-RR(lc+1)*RR(lc));
  }
  /*
   * bottom layer
   */
  irow++;
  lc = ncut;
  /*
   * SUPERD(irow) = undefined
   */
  wk = exp(-KK(lc)*(TAUCPR(lc)-TAUCPR(lc-1)));
  if (lyrcut) {
    SUBD(irow) = 1.;
    DIAG(irow) = RR(lc)*wk;
  }
  else {
    SUBD(irow) = 1.-2.*ds->bc.albedo*cmu*RR(lc);
    DIAG(irow) = (RR(lc)-2.*ds->bc.albedo*cmu)*wk;
  }

  /*
   * NOTE: If not allocating memory with swappablecalloc(), need to zero out b here.
   */

  /*
   * Construct -b-, for parallel beam + bottom reflection + thermal emission at top and/or bottom
   *
   * Top boundary, right-hand-side of eq. KST(28)
   */
  lc   = 1;
  irow = 1;
  B(irow) = -YB_0D(lc)-YP_0D(lc)+ds->bc.fisot+tplanck;
  /*
   * Continuity condition for layer interfaces, right-hand-side of eq. KST(29)
   */
  for (lc = 1; lc <= nloop; lc++) {
    fact1     = exp(-ZB_A(lc+1)*TAUCPR(lc));
    fact2     = exp(-ZP_A(lc+1)*TAUCPR(lc));
    fact3     = exp(-ZB_A(lc  )*TAUCPR(lc));
    fact4     = exp(-ZP_A(lc  )*TAUCPR(lc));
    rpp1_m    = fact1*(YB_0D(lc+1)+YB_1D(lc+1)*TAUCPR(lc))+fact2*(YP_0D(lc+1)+YP_1D(lc+1)*TAUCPR(lc));
    rp_m      = fact3*(YB_0D(lc  )+YB_1D(lc  )*TAUCPR(lc))+fact4*(YP_0D(lc  )+YP_1D(lc  )*TAUCPR(lc));
    rpp1_p    = fact1*(YB_0U(lc+1)+YB_1U(lc+1)*TAUCPR(lc))+fact2*(YP_0U(lc+1)+YP_1U(lc+1)*TAUCPR(lc));
    rp_p      = fact3*(YB_0U(lc  )+YB_1U(lc  )*TAUCPR(lc))+fact4*(YP_0U(lc  )+YP_1U(lc  )*TAUCPR(lc));
    B(++irow) = rpp1_p-rp_p-RR(lc+1)*(rpp1_m-rp_m);
    B(++irow) = rpp1_m-rp_m-RR(lc  )*(rpp1_p-rp_p);
  }
  /*
   * Bottom boundary
   */
  lc = ncut;
  if (lyrcut) {
    /*
     * Right-hand-side of eq. KST(30)
     */
    B(++irow) = -exp(-ZB_A(ncut)*TAUCPR(ncut))*(YB_0U(ncut)+YB_1U(ncut)*TAUCPR(ncut))
                -exp(-ZP_A(ncut)*TAUCPR(ncut))*(YP_0U(ncut)+YP_1U(ncut)*TAUCPR(ncut));
  }
  else {
    sum = cmu*ds->bc.albedo*(exp(-ZB_A(ncut)*TAUCPR(ncut))*(YB_0D(ncut)+YB_1D(ncut)*TAUCPR(ncut))
                            +exp(-ZP_A(ncut)*TAUCPR(ncut))*(YP_0D(ncut)+YP_1D(ncut)*TAUCPR(ncut)));
   if (ds->bc.umu0 <= 0.) {
     refflx = 0.;
   }
   else {
     refflx = 1.;
   }
   B(++irow) = 2.*sum+ds->bc.albedo*ds->bc.umu0*ds->bc.fbeam/M_PI*refflx*EXPBEA(ncut)+(1.-ds->bc.albedo)*bplanck
               -exp(-ZB_A(ncut)*TAUCPR(ncut))*(YB_0U(ncut)+YB_1U(ncut)*TAUCPR(ncut))
               -exp(-ZP_A(ncut)*TAUCPR(ncut))*(YP_0U(ncut)+YP_1U(ncut)*TAUCPR(ncut));

 }
 /*
  * solve for constants of integration by inverting matrix KST(38-41)
  */
  nrow = irow;

  /*
   * NOTE: If not allocating memory with swappablecalloc(), need to zero out cband here.
   */

  for (irow = 1; irow <= nrow; irow++) {
    CBAND(1,irow) = 0.;
    CBAND(3,irow) = DIAG(irow);
  }
  for (irow = 1; irow <= nrow-1; irow++) {
    CBAND(2,irow+1) = SUPERD(irow);
  }
  for (irow = 2; irow <= nrow; irow++) {
    CBAND(4,irow-1) = SUBD(irow);
  }

  c_sgbfa(cband,(9*(ds->nstr/2)-2),nrow,1,1,ipvt,&info);
  job = 0;
  c_sgbsl(cband,(9*(ds->nstr/2)-2),nrow,1,1,ipvt,b,job);

  /*
   * unpack
   */
  irow = 0;
  for (lc = 1; lc <= ncut; lc++) {
    /* downward direction */
    LL(1,lc) = B(++irow);

    /* upward direction */
    LL(2,lc) = B(++irow);
  }

  return;
}

/*============================= end of c_twostr_solve_bc() ===============*/
