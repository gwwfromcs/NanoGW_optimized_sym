!
! Weiwei Gao, Feb. 2018
!
! In the density fitting method, the products of wave functions are
! written as linear combinations of some interpolation vectors zeta(r).
! See equation 5 in J. Chem. Theory. Comput. 2017, 13, 5420-5431:
!
! \phi_i(r)\psi_j(r) = \sum^{Nu}_{u=1} \zeta_u(r) \phi_i(r_u) \psi_j(r_u)
! 
! <==> written in a matrix form: Z = Zeta * C
! __________________________________________________________________________________
! |  Matrix      |     Shape      |        Contents                                |
! |______________|________________|________________________________________________|
! |   P(r,r_u)   |  Ng * N_u      |    \sum_i \phi_i(r)\psi_i(r_u)                 |
! |              |                |     i=1...Nv, r are full-grid points           | 
! |______________|________________|________________________________________________|
! |   Q(r,r_u)   |  Ng * N_u      |    \sum_i \phi_i(r)\psi_i(r_u)                 |
! |              |                |     i=1...Nc, r are full-grid points           | 
! |______________|________________|________________________________________________|
! |  Zeta_u(r)   |  Ng *  Nu      |    {..., \zeta_u(r), ...}                      |
! |              |                |     u=1...Nu                                   |
! |______________|________________|________________________________________________|
! |  C(r_u,i,j)  |  Nu * (Nc*Nv)  |    {..., \phi_i(ru)*\psi_j(ru), ...}           |
! |              |                |                                                |
! |______________|________________|________________________________________________|
! _________________________________________________________________________
! |  Matrix      |                  Parallelization                       |
! |______________|________________________________________________________|
! |   P(r,r_u)   |             Each proc store part of Z                  |
! |              |             P(ngfl, Nc*Nv), ngfl = mydim               |
! |______________|________________________________________________________|
! |   Q(r,r_u)   |             Each proc store part of Z                  |
! |              |             Q(ngfl, Nc*Nv), ngfl = mydim               |
! |______________|________________________________________________________|
! |  Zeta_u(r)   |             Each proc store part of zeta               |
! |              |         zeta(ngfl, Nu), ngfl = mydim*ntrans            |
! |______________|________________________________________________________|
! |  C(r_u,i,j)  |             Every proc store a full copy of C          |
! |              |              Cmtrx ( Nu, Nc*Nv )                       |
! |______________|________________________________________________________|
!
!
! This subroutine calculates the interpolation vectors zeta_u(r)
! 
! n_intp   : the number of interpolation vectors or points, or Nu
! n_intp_r : the number of interpolation vectors or points in reduced r-space domain
! intp     : the index of interpolation points in the full grid
! zeta(gvec%nr, n_intp_r) : the interpolation vectors 
! kpt%wfn(isp,ikp)%dwf(:,:) : store the wavefunctions \phi_i, \psi_j
!
! For now, this is only partially parallelized
! 1. each processor evaluates part of the matrix Z and part of the matrix C
! 2. collect matrix C from different procs. distribute matrix Z to procs.
! 3. every proc in w_grp solve part of the linear equations to get zeta
!
! This subroutine consists of the following main steps:
!
! Step 1.0: Prepare P(r,r_u)
! Step 1.5: Prepare A=C.C^T=P(r_u,r_v)Q(r_u,r_v), B=Z.C^T=P(r,r_u)Q(r,r_v)
! Step 1.6: Deallocate P and Q
! Step 2.0: solve a linear quation A.zeta=B to get zeta
! Step 3.0: Calculate Mmtrx(1) = <zeta(i)| V_coul(r,r') | zeta(j)> and
!       Mmtrx(2) = <zeta(i)| f^LDA(r) | zeta(j)>
! Step 3.5: Calculate C(r_u,i,j)

subroutine isdf_parallel_sym( gvec, pol_in, kpt, n_intp_r, intp_r, &
      nspin, ncv, maxncv, invpairmap, nv, maxnv, ivlist, nc, maxnc, iclist, kflag, &
      Cmtrx, Mmtrx, verbose )
 
  use typedefs
  use mpi_module
  use myconstants
  use fft_module
  use xc_functionals
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  type(gspace), intent(in) :: gvec
  type(polinfo), intent(in), dimension(2) :: pol_in
  type(kptinfo), intent(in) :: kpt
  ! n_intp_r: number of interpolation points in reduced r-space domain
  ! intp_r(n_intp_r): the index of interpolation points in full grids
  ! ncv(2): the number of wave function pairs for two spins
  !
  ! invpairmap maps the index of pair |cv> to c and v
  integer, intent(in) :: n_intp_r, intp_r(n_intp_r),          &
    nspin,                                                    &
    ncv(nspin,kpt%nk,gvec%syms%ntrans),                       &
    maxncv,                                                   &
    invpairmap(2, maxncv, nspin, kpt%nk, gvec%syms%ntrans),   &
    nv(nspin, kpt%nk, gvec%syms%ntrans),                      &
    maxnv,                                                    &
    ivlist(maxnv, nspin, kpt%nk, gvec%syms%ntrans),           & 
    nc(nspin, kpt%nk, gvec%syms%ntrans),                      &
    maxnc,                                                    &
    iclist(maxnc, nspin, kpt%nk, gvec%syms%ntrans), kflag
  
  ! If we consider periodic systems, we should change these from real(dp) to SCALAR
  real(dp), intent(out) :: Cmtrx(n_intp_r, maxncv,   nspin, kpt%nk, gvec%syms%ntrans)  
  real(dp), intent(out) :: Mmtrx(n_intp_r, n_intp_r, nspin, nspin, kpt%nk, 2, gvec%syms%ntrans)

  ! If verbose is true, then print out additional debug information
  logical, intent(in) :: verbose

  ! ------------------ Local variables ------------------
  ! P(gvec%mydim, n_intp_r, nspin, kpt%nk)
  ! Q(gvec%mydim, n_intp_r, nspin, kpt%nk)
  ! Cmtrx(n_intp_r, Nc*Nv, nspin, kpt%nk)
  ! Issue need to be addressed: if there are N wfn_grp, does it mean these matrix need to be
  ! calculated by all processors at the same time, or can we distribute the workload
  ! later ??
  real(dp), allocatable ::            &
     PsiV(:,:,:), PsiV_intp(:,:,:),   &   ! PsiV: wfn on reduced domain
     PsiC(:,:,:), PsiC_intp(:,:,:),   &
     P(:,:,:), P_intp(:,:,:),         &   ! P on reduced domain 
     Q(:,:,:), Q_intp(:,:,:),         &   ! Q on reduced domain 
     zeta(:,:,:,:), tmp_Zmtrx(:), Zmtrx(:), Zmtrx_loc(:), &
     fxc(:,:,:), fxc_loc(:,:,:), fzeta(:), &
     ! matrices and vectors used for solving linear equations
     Amtrx(:,:), Bmtrx(:,:), Xmtrx(:,:), tmpmtrx(:,:,:,:), &
     rho_h(:), rho_h_distr(:)
  real(dp) :: diff, weight, qkt(3), tsec(2), &
     vcharac(gvec%syms%ntrans), ccharac(gvec%syms%ntrans)
  integer, allocatable :: inv_ivlist(:,:,:,:), inv_iclist(:,:,:,:)
  ! counters and temporary integers
  integer :: ipt, ii, jj, iv, ic, icv, irp, jrp, rsp, csp, i_row, i_col, lrp1, lrp2, & 
             IVV, ICC, JVV, JCC, isp, ikp, errinfo, ipe, &
             ivrp, icrp, maxivv, maxicc
  integer :: status, mpi_status(MPI_STATUS_SIZE)
  ! The index of w_grp, r_grp, processor that works on the calculation of
  ! zeta(:,:,nspin,nk)
  integer :: workwgrp, workrgrp, workproc
  ! Each processor in a w_grp store part of the wave function
  ! offset: index of grid point from which a processor start to store the wave function
  ! ncount: number of elements of wave functions that a processor stores
  integer, dimension(0:w_grp%npes-1) :: offset, ncount
  ! temporary dummy variable
  integer, dimension(0:w_grp%npes-1) :: idum
  ! the number of grid points in irreducible wedge, ngr = gvec%nr
  integer :: ngr, ngrid
  ! the number of full grid point, ngf = ngr * (# of sym operations)
  integer :: ngfl, iptf, iptr, ioff, ioff1, ioff2

  ! variables for debug and test of accuracy 
  character(50) :: dbg_filename = "isdf_dbg.dat"
  integer :: dbgunit = 20171130
  ! external functions
  real(dp), external :: ddot
  type(xc_type) :: xc_lda
  !
  workrgrp = 0
  workwgrp = 0
  call MPI_BCAST( workrgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo )
  call MPI_BCAST( workwgrp, 1, MPI_INTEGER, peinf%masterid, peinf%comm, errinfo )
  !
  if( w_grp%master .and. verbose ) then
     !
     write(*,*) "call isdf_parallel(), write debug info to ", dbg_filename
     !
     open(dbgunit, file=dbg_filename, form='formatted', status='replace')
     !
  endif
  !
  ngrid= w_grp%mydim
  ! the number of real-space grid points in reduced real-space domain
  ngr = gvec%nr 
  ! the number of grid points of interpolation vectors zeta(:) stored in each processor
  ngfl = w_grp%mydim * gvec%syms%ntrans
  !
  ! each processor stores part of P, Q and zeta
  !
  ! Allocating intermediate variables
  ALLOCATE(P     (w_grp%mydim, n_intp_r, gvec%syms%ntrans ))
  ALLOCATE(Q     (w_grp%mydim, n_intp_r, gvec%syms%ntrans ))
  ALLOCATE(P_intp(n_intp_r,    n_intp_r, gvec%syms%ntrans ))
  ALLOCATE(Q_intp(n_intp_r,    n_intp_r, gvec%syms%ntrans ))
  ALLOCATE(Amtrx (n_intp_r,    n_intp_r                   ))
  ALLOCATE(Bmtrx (n_intp_r,    w_grp%mydim                ))
  ALLOCATE(Xmtrx (n_intp_r,    w_grp%mydim                ))
  ALLOCATE(zeta  (w_grp%mydim, n_intp_r, nspin, gvec%syms%ntrans )) ! interpolation vectors
  ALLOCATE(PsiV( w_grp%mydim,   maxnv, gvec%syms%ntrans))
  ALLOCATE(PsiV_intp( n_intp_r, maxnv, gvec%syms%ntrans ))
  ALLOCATE(PsiC( w_grp%mydim,   maxnc, gvec%syms%ntrans ))
  ALLOCATE(PsiC_intp( n_intp_r, maxnc, gvec%syms%ntrans ))
  !
  ! initialize matrices with zero
  !
  Cmtrx = zero
  zeta  = zero
  !
  if ( w_grp%master .and. verbose ) then
     write(dbgunit, '(a)') " Index of interpolation points in full domain: "
     write(dbgunit, '(5i15)') ( intp_r(ii), ii=1, n_intp_r )
  endif
  !
  idum = 0 
  idum(w_grp%inode) = w_grp%offset + 1
  call MPI_ALLREDUCE(idum, offset, w_grp%npes, MPI_INTEGER, MPI_SUM, &
     w_grp%comm, errinfo)
  !
  idum = 0
  idum(w_grp%inode) = w_grp%mydim
  call MPI_ALLREDUCE(idum, ncount, w_grp%npes, MPI_INTEGER, MPI_SUM, &
     w_grp%comm, errinfo)
  !
  if ( verbose .and. w_grp%master ) then
     !
     write(dbgunit, *) " in isdf() "
     write(dbgunit, *) " w_grp%mygr ", w_grp%mygr, " w_grp%inode = ", w_grp%inode, &
     " workwgrp = ", workwgrp, &
     " workrgrp = ", workrgrp, &
     " offset: ", ( offset(ii), ii=0, w_grp%npes-1), &
     " ncount: ", ( ncount(ii), ii=0, w_grp%npes-1)
     !
  endif
  !
  maxivv = maxval(ivlist)
  maxicc = maxval(iclist)
  ! ivlist(:) maps the index of valence states used in 
  !   calculation (i.e., stored in memory) to the real index of valence states
  allocate(inv_ivlist(maxivv, nspin, kpt%nk, gvec%syms%ntrans)) ! For now, only deal with confined system. So we assume kpt%nk=1, and there is no dependence on k here
  allocate(inv_iclist(maxicc, nspin, kpt%nk, gvec%syms%ntrans))
  ! inv_ivlist(:) maps the real index of valence states 
  !   to the index of valence states used in the calculation
  inv_ivlist = 0
  inv_iclist = 0
  !
  ! First, you have to select the interpolation points in the reduced domain from intp(:)
  do ikp = 1, kpt%nk
     do isp = 1, nspin
        !
        do irp = 1, gvec%syms%ntrans
           do iv = 1, nv(isp,ikp,irp)
              IVV = ivlist(iv, isp, ikp, irp)
              inv_ivlist(IVV, isp, ikp, irp) = iv
           enddo
           !
           do ic = 1, nc(isp,ikp,irp)
              ICC = iclist(ic, isp, ikp, irp)
              inv_iclist(ICC, isp, ikp, irp) = ic
           enddo
           !
        enddo
        !
     enddo
  enddo 
  !
  Mmtrx = zero ! Mmtrx dimension: Mmtrx(n_intp_r, n_intp_r, nspin, nspin)
  !
  ! qkt is set to zero. This is only valid for
  !  tests of nonperiodic system, and should be updated later.
  !
  qkt = 0
  ALLOCATE(rho_h(ngr))     ! note: gvec%nr is equal to w_grp%nr
  ALLOCATE(rho_h_distr(w_grp%ldn*w_grp%npes))
  !
  ! kflag = 0 : calculate kernel K^x  ( Coulomb )
  !         1 :                  K^x + K^f ( Coulomb + 1/2 * F_xc )
  !         2 :                  K^f   ( 1/2 * F_xc )
  !         3 :                  K^t   ( 1/2 * F_xc for spin triplet )
  !
  ! Generate Coulomb potential
  !
  if (kflag < 2 ) then
     call dinitialize_FFT(peinf%inode, fft_box)
     call dcreate_coul_0D(gvec%bdot, qkt, fft_box)
  endif
  !
  ! Calculate LDA kernel
  !
  if ( kflag > 0 ) then
     !
     ALLOCATE( fxc( w_grp%mydim, nspin, nspin ), stat = errinfo )
     ALLOCATE( fzeta(ngfl), stat = errinfo )
     fxc = zero
     ! 
     ! Copy the charge density to fxc
     !
     do isp = 1, nspin
        call dcopy( w_grp%mydim, kpt%rho(w_grp%offset+1, isp), 1, &
           fxc(1, isp, 1), 1 )
     enddo 
     call xc_init( nspin, XC_LDA_X, XC_LDA_C_PZ, 0, zero, one, .false., xc_lda )
     xc_lda%has_grad = .false.
     call fxc_get( xc_lda, nspin, w_grp%mydim, kflag, fxc )
     call xc_end( xc_lda )
     !
  endif
  !
  do ikp = 1, kpt%nk
     !
     do isp = 1, nspin
        !
        ! The following loop calculate P(r,r_u,jrp), Q(r,r_u,jrp) for all representations 
        !
        do jrp = 1, gvec%syms%ntrans
           !
           ! initialize matrices with zero
           !
           P = zero 
           Q = zero
           PsiV  = zero
           PsiC  = zero
           PsiV_intp = zero
           PsiC_intp = zero
           !
           if (verbose .and. w_grp%master) &
              write(dbgunit, *) ' isp = ', isp, ', ikp = ', ikp
           !
           do iv = 1, nv(isp,ikp,jrp)
              IVV = ivlist(iv,isp,ikp,jrp)
              JVV = kpt%wfn(isp,ikp)%map(IVV)
              ! PsiV_i(r)
              call dcopy(w_grp%mydim, & 
                kpt%wfn(isp,ikp)%dwf(1,JVV),1,PsiV(1,iv,jrp),1)
           enddo
           !
           do ic = 1, nc(isp,ikp,jrp)
              ICC = iclist(ic,isp,ikp,jrp)
              JCC = kpt%wfn(isp,ikp)%map(ICC)
              ! PsiC_i(r)
              call dcopy(w_grp%mydim, & 
                kpt%wfn(isp,ikp)%dwf(1,JCC),1,PsiC(1,ic,jrp),1)
           enddo
           !
           ! pick the interpolation points
           !
           do ipt = 1, n_intp_r
              iptf = intp_r(ipt)
              iptr = iptf / gvec%syms%ntrans
              ioff1 = offset(w_grp%inode)
              ioff2 = offset(w_grp%inode+1)
              if ( iptr .ge. ioff1 &
                   .and. &
                   iptr .lt. ioff2 ) then
                jj = iptr - ioff1 + 1
                ! if ( icv .eq. 1 .and. verbose ) write(*, '(4(a,i5))') &
                !    ", wgrp%inode ", w_grp%inode, ", ipt ", ipt, &
                !    ", iptf ", iptf, ", jj ", jj
                PsiV_intp(ipt, 1:nv(isp,ikp,jrp), jrp) = PsiV(jj, 1:nv(isp,ikp,jrp), jrp)
                PsiC_intp(ipt, 1:nc(isp,ikp,jrp), jrp) = PsiC(jj, 1:nc(isp,ikp,jrp), jrp)
              endif
           enddo
           !
           call MPI_ALLREDUCE(MPI_IN_PLACE, PsiV_intp(1,1,jrp), n_intp_r*nv(isp,ikp,jrp), MPI_DOUBLE, MPI_SUM, &
             w_grp%comm, errinfo)
           call MPI_ALLREDUCE(MPI_IN_PLACE, PsiC_intp(1,1,jrp), n_intp_r*nc(isp,ikp,jrp), MPI_DOUBLE, MPI_SUM, &
             w_grp%comm, errinfo)
           !
           ! Prepare P and Q
           !
           ! P(r,r_u,jrp) = \sum_{V~jrp} \PsiV(r) \PsiV(r_u)
           call dgemm('n','t',w_grp%mydim, n_intp_r, nv(isp,ikp,jrp), one, &
             PsiV(1,1,jrp), w_grp%mydim, PsiV_intp(1,1,jrp), n_intp_r, zero, P(1,1,jrp), w_grp%mydim)
           ! Q(r,r_u,jrp) = \sum_{C~jrp} \PsiC(r) \PsiC(r_u)
           call dgemm('n','t',w_grp%mydim, n_intp_r, nc(isp,ikp,jrp), one, &
             PsiC(1,1,jrp), w_grp%mydim, PsiC_intp(1,1,jrp), n_intp_r, zero, Q(1,1,jrp), w_grp%mydim)
           ! Calculate P(r_u,r_u,jrp), and Q(r_u,r_u,jrp)
           call dgemm('n','t',n_intp_r,n_intp_r, nv(isp,ikp,jrp), one, &
             PsiV_intp(1,1,jrp), n_intp_r, PsiV_intp(1,1,jrp), n_intp_r, zero, &
             P_intp(1,1,jrp), n_intp_r)
           call dgemm('n','t',n_intp_r,n_intp_r, nc(isp,ikp,jrp), one, &
             PsiC_intp(1,1,jrp), n_intp_r, PsiC_intp(1,1,jrp), n_intp_r, zero, &
             Q_intp(1,1,jrp), n_intp_r)

           !
        enddo ! jrp loop
        do jrp = 1, gvec%syms%ntrans
           !
           ! Prepare Cmtrx
           !
           do icv = 1, ncv(isp,ikp,jrp)
              IVV  = invpairmap(1,icv,isp,ikp,jrp)
              ivrp = kpt%wfn(isp,ikp)%irep(IVV)
              ICC  = invpairmap(2,icv,isp,ikp,jrp)
              icrp = kpt%wfn(isp,ikp)%irep(ICC)
              Cmtrx(1:n_intp_r, icv, isp, ikp, jrp) = PsiV_intp(1:n_intp_r, inv_ivlist(IVV,isp,ikp,ivrp), ivrp) * &
                PsiC_intp(1:n_intp_r, inv_iclist(ICC,isp,ikp,icrp), icrp) ! This is element-wise multiplication
           enddo
        enddo
        !
        ! Calculate zeta(r,n_intp_r,jrp) for all representations
        !
        do jrp = 1, gvec%syms%ntrans/r_grp%num
           irp = r_grp%g_rep(jrp)
           Amtrx = zero
           Bmtrx = zero
           Xmtrx = zero
           ! Calculate A and B
           ! Loop over all the reps.
           do lrp1 = 1, gvec%syms%ntrans
              !
              ! The direct product of lrp1 and lrp2 should be irp: lrp1 * lrp2 = irp
              !
              do lrp2 = 1, gvec%syms%ntrans
                 if(gvec%syms%prod(lrp1,lrp2) == irp) exit
              enddo
              !
              ! ---------------------
              ! For each irp, spin and ikp, calculate zeta
              ! zeta dimension: Ng*Nu
              ! Set A = ( sum_{lrp1} P_intp(r_u,r_u,lrp1).Q_intp(r_u,r_u,lrp2) )^T  dimension: n_intp_r * n_intp_r
              !     B = ( sum_{lrp1} P(r,r_u,lrp1).Q(r,r_u,lrp2) )^T                dimension: n_intp_r * w_grp%mydim
              ! A is a symmetric matrix!
              ! we solve the linear equation A*(zeta^T) = B to get zeta^T for representation irp
              ! zeta^T = A^-1 * B = (C * C^T)^-1 * C * Z^T
              !
              !  Matrix dimensions:
              !  zeta(ngfl, n_intp_r, :, :)
              !  Amtrx(n_intp_r, n_intp_r)       intermediate variable, store C * C^T
              !  Bmtrx(n_intp_r, w_grp%mydim)  intermediate variable, store C * Z^T
              !  Xmtrx(n_intp_r, w_grp%mydim)         intermediate variable, store zeta^T
              ! ---------------------
              !
              ! calculate A = P_intp.Q_intp (Note: This is an element-wise multipliation)
              !
              Amtrx(1:n_intp_r,1:n_intp_r) = P_intp(1:n_intp_r,1:n_intp_r,lrp1) * &
                                             Q_intp(1:n_intp_r,1:n_intp_r,lrp2)
              !
              if ( verbose .and. w_grp%master ) then
                 write(dbgunit,'(" A = C*C^T = P_intp*Q_intp = ")')
                 call printmatrix ( Amtrx(1,1), n_intp_r, n_intp_r, dbgunit )
              endif
              !
              ! calculate B = (P.Q)^T (Note: This is an element-wise multiplication)
              !
              Bmtrx(1:n_intp_r,1:w_grp%mydim) = transpose( P(1:w_grp%mydim, 1:n_intp_r, lrp1) * &
                                              Q(1:w_grp%mydim, 1:n_intp_r, lrp2) )
              !
              if ( verbose .and. w_grp%master ) then
                 write(dbgunit,'(" B = C*Z^T = P*Q = ")')
                 call printmatrix ( Bmtrx(1,1), n_intp_r, w_grp%mydim, dbgunit )
              endif
              !
           enddo ! lrp1
           !
           ! solver the linear equation A * X = B
           !
           call dlinear_solver( n_intp_r, w_grp%mydim, Amtrx, Bmtrx, Xmtrx, w_grp%inode, verbose )
           !
           ! Copy Xmtrx to zeta
           !
           ! if ( verbose ) write(*,*) 'zeta', isp,ikp
           do ii = 1, n_intp_r
              do jj = 1, w_grp%mydim
                 zeta( jj, ii, isp, irp ) = Xmtrx( ii, jj )
              enddo ! jj loop
           enddo ! ii loop
           !
        enddo ! jrp loop
        !
     enddo ! isp loop
     !
     do jrp = 1, gvec%syms%ntrans/r_grp%num
        ! jrp is the index of representation belong to this r_grp
        ! irp is the real index of the representation
        irp = r_grp%g_rep(jrp)
        !
        ! Now calculate <zeta_u(r,ispin)|V(r,r')|zeta_w(r',jspin)>, where u, w = 1, ..., n_intp
        !  and store it in Mmtrx(n_intp_r, n_intp_r)
        !
        do rsp = 1, nspin
           do csp = 1, nspin
              ! 
              ! Exp: n_intp_r=1000, r_grp%npes=10, w_grp%npes=5
              !      then w_grp%mygr=0 work on ii = 1,2,3,4,5, 11,12,13,14,15, 21,22,23,24,25, ...
              !           w_grp%mygr=1 work on ii = 5,6,7,8,9, 16,17,18,19,20, 26,27,28,29,30, ...
              do i_row = 1, n_intp_r, r_grp%npes
                 !
                 ! initialize rho_h and rho_h_local
                 !
                 rho_h = 0.d0
                 rho_h_distr = 0.d0
                 !
                 do ipe = 0, w_grp%npes-1
                    ii   = w_grp%mygr*w_grp%npes+i_row+ipe
                    ! Note: w_grp%ldn = w_grp%mydim or w_grp%mydim+1
                    ioff = w_grp%ldn*ipe + 1
                    if (ii > n_intp_r) cycle
                    ! Note:  ngrid = w_grp%mydim
                    ! rho_h_distr are distributed over w_grp
                    call dcopy( ngrid, &
                      zeta( 1, ii, rsp, irp ), 1, &
                      rho_h_distr( ioff ), 1 )
                    if ( kflag > 0 ) then
                       call dcopy( ngrid, zeta( 1, ii, rsp, irp ), 1, &
                         fzeta(ioff), 1 )
                       ! calculate fzeta = f_lda * zeta
                       call dmultiply_vec( ngrid, fxc(1,rsp,csp), fzeta(ioff) ) 
                    endif ! kflag > 0
                 enddo ! ipe loop
                 !
                 ii = w_grp%mygr*w_grp%npes + i_row + w_grp%inode
                 if ( kflag < 2 ) then
                    if ( ii < n_intp_r) then
                    call dgather(1,rho_h_distr,rho_h)
                    !
                    ! solve poisson equation to get: rho_h(r) = \int V_c(r,r') zeta_ii(r') dr'
                    !
                    call dpoisson(gvec, rho_h, irp)
                    call dscatter(1,rho_h_distr,rho_h)
                    endif
                 endif ! kflag < 2 
                 !
                 ! Mmtrx is a symmetric matrix
                 !
                 do ipe = 0, w_grp%npes-1
                    ! ii is real row index
                    ii = w_grp%mygr*w_grp%npes+i_row+ipe
                    ioff=w_grp%ldn*ipe+1
                    do i_col = 1, n_intp_r
                      !
                      if (rsp == csp .and. i_col < ii) cycle
                      !
                      ! calculate: \int zeta_jj(r) V_c(r,r') zeta_ii(r') drdr'
                      !          = \int zeta_jj(r) rho_h(r) dr
                      !
                      if ( kflag > 0 ) then
                        Mmtrx( ii, i_col, rsp, csp, ikp, 2, jrp ) = &
                          ddot( ngrid, zeta( 1, i_col, csp, jrp ), 1, fzeta(ioff), 1 ) / gvec%hcub
                        !
                        ! Mmtrx is a symmetric matrix
                        !
                        Mmtrx( i_col, ii, csp, rsp, ikp, 2, jrp ) = &
                          Mmtrx( ii, i_col, rsp, csp, ikp, 2, jrp )
                      endif
                      if ( kflag < 2 ) then
                        Mmtrx( ii, i_col, rsp, csp, ikp, 1, jrp ) = &
                          ddot( ngrid, zeta( 1, i_col, csp, jrp ), 1, rho_h_distr(ioff), 1 ) / gvec%hcub
                        !
                        ! Mmtrx is a symmetric matrix
                        !
                        Mmtrx( i_col, ii, csp, rsp, ikp, 1, jrp) = &
                          Mmtrx( ii, i_col, rsp, csp, ikp, 1, jrp )
                      endif
                      !
                    enddo ! i_col loop
                 enddo ! ipe loop
                 !
              enddo ! i_row loop 
              !
           enddo ! csp loop
           !
        enddo ! rsp loop
        !
     enddo ! jrp loop
     !
  enddo ! ikp loop
  DEALLOCATE( PsiV )
  DEALLOCATE( PsiV_intp )
  DEALLOCATE( PsiC )
  DEALLOCATE( PsiC_intp )
  !
  ! clean up all the allocated variables
  !
  if ( peinf%master .and. verbose ) then
     write(*,*) " DEALLOCATING arrays"
  endif
  DEALLOCATE( Amtrx )
  DEALLOCATE( Bmtrx )
  DEALLOCATE( Xmtrx )
  DEALLOCATE( P )
  DEALLOCATE( P_intp )
  DEALLOCATE( Q )
  DEALLOCATE( Q_intp )
  DEALLOCATE( inv_ivlist )
  DEALLOCATE( inv_iclist )
  !
  if ( kflag > 0 ) then
    DEALLOCATE(fxc)
    DEALLOCATE(fzeta)
  endif
  DEALLOCATE(zeta) ! no longer needed
  DEALLOCATE(rho_h)
  DEALLOCATE(rho_h_distr)
  !
  if ( kflag > 0 ) then
     do jrp = 1, gvec%syms%ntrans
        call MPI_ALLREDUCE( MPI_IN_PLACE, Mmtrx(1,1,1,1,1,2,jrp), n_intp_r * n_intp_r * nspin * nspin * kpt%nk, &
          MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
     enddo
  endif
  if ( kflag < 2 ) then
     do jrp = 1, gvec%syms%ntrans
        call MPI_ALLREDUCE( MPI_IN_PLACE, Mmtrx(1,1,1,1,1,1,jrp), n_intp_r * n_intp_r * nspin * nspin * kpt%nk, &
          MPI_DOUBLE, MPI_SUM, w_grp%comm, errinfo )
     enddo
  endif
  !
  if ( workrgrp == r_grp%mygr .and. workwgrp == w_grp%mygr ) then
     if ( w_grp%master .and. verbose ) then
        write( *, * ) " ikp = ", ikp, " Mmtrx (:, :, isp=1, isp=1, ikp=1, 1) = "
        call printmatrix ( Mmtrx (1,1,1,1,1,1,jrp), n_intp_r, n_intp_r, dbgunit ) 
        close ( dbgunit )
     endif
  endif
  return
  !
end subroutine isdf_parallel_sym
