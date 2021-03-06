!===================================================================
!
! Diagonalizes the TDLDA eigenvalue equations.
!
! Input file is "tdlda.in" (optional). Output
! is writen on screen and on files "eigenvalues_rpa" (RPA-LDA energy
! eigenvalues and oscillator strengths), "eigenvalues_lda" (the same
! for TDLDA) and "pol_eig.dat" (eigenvectors).
!
! This code follows the numerical methodology presented in:
!    M.E. Casida, in "Recent Advances in Density-Functional Methods", ed. D.P. Chong (1995)
!    I. Vasiliev, S. Ogut and J.R. Chelikowsky, Phys. Rev. B 65, 115416 (2002)
!
! Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
! mtiago@ices.utexas.edu
!
! First version written by Murilo Tiago, Univ. of Minnesota, April 2004
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 1, or (at your option)
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA
!
!-------------------------------------------------------------------
program tdlda
  use typedefs
  use mpi_module
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  type (gspace) gvec
  type (kptinfo) :: kpt
  type (qptinfo) :: qpt
  type (polinfo), dimension(2) :: pol_in
  type (polinfo), allocatable :: pol(:,:)
  type (kernelinfo), allocatable :: k_p(:,:)

  character (len=800) :: lastwords
  character (len=40), allocatable :: routnam(:)
  logical :: nolda, tamm_d, rpaonly, trip_flag, noxchange, trunc_c, init_gr
  integer :: ii, jj, isp, irp, iq, nspin, nmap, nbuff, lcache, dft_code
  real(dp) :: tsec(2), mem1, xsum, xmax, tdldacut
  logical, allocatable :: wmap(:)
#ifdef MPI
  integer :: info
#endif

  ! W Gao dbg
  logical :: doisdf
  real(dp), allocatable :: Cmtrx(:,:,:,:,:), Mmtrx(:,:,:,:,:,:,:)
  ! the number of val and cond states pairs. Since there could be two spins,
  ! ncv(1) correspond to spin up, and ncv(2) correspond to spin down.
  ! Assumption: for different kpts, the number of states with the same 
  !  spin are the same 
  integer :: n_intp, n_intp_r, maxicc, maxivv, iv, ic, &
           ivv, icc, maxnc, maxnv, maxncv, &
           ihomo, ikp, intp_type, isdf_type, &
           kflag, vcirp, virp, cirp, ipt
  ! cvt.f90
  integer, allocatable :: intp(:), pairmap(:,:,:,:,:), invpairmap(:,:,:,:,:), icv(:), &
           ncv(:,:,:), ivlist(:,:,:,:), iclist(:,:,:,:), nv(:,:,:), nc(:,:,:), timerlist(:), &
           intp_r(:)

  ! WG debug
  integer :: outdbg
  logical :: verbose
  character (len=20) :: dbg_filename

  !-------------------------------------------------------------------
  ! Initialization.
  !
  call header('TDLDA')
  !-------------------------------------------------------------------
  ! Read input info.
  !
  if(peinf%master) write (*,*) " Call input_g"
  call input_g(pol_in,qpt,tdldacut,nbuff,lcache,w_grp%npes,&
       nolda,tamm_d,r_grp%num,dft_code,doisdf,n_intp,intp_type,isdf_type,.false.)
  if(peinf%master) write (*,*) " input_g done"
  call input_t(tamm_d,rpaonly,trip_flag,noxchange,trunc_c)
  
  ! W Gao open dbg files
  write(dbg_filename,"(i7)") peinf%inode
  outdbg = peinf%inode+198812
  dbg_filename = "kernel_dbg"//adjustl(dbg_filename)
  open(outdbg,file=dbg_filename,status='unknown',iostat=info) 

  !-------------------------------------------------------------------
  ! Determine the set of wavefunctions to read: if n-th wavefunction is
  ! needed, then wmap(n) = .true.; otherwise, wmap(n) = .false.
  !
  nmap = 0
  ! if any of pol_in(1)%ncond, .... is larger than 0, then
  if (max(pol_in(1)%ncond,pol_in(1)%nval, &
       pol_in(2)%ncond,pol_in(2)%nval) > 0) then
     ! nmap is equal to the max of pol_in(1)%cmap(:), ...
     if (pol_in(1)%ncond > 0) nmap = max(nmap,maxval(pol_in(1)%cmap))
     if (pol_in(1)%nval > 0) nmap = max(nmap,maxval(pol_in(1)%vmap))
     if (pol_in(2)%ncond > 0) nmap = max(nmap,maxval(pol_in(2)%cmap))
     if (pol_in(2)%nval > 0) nmap = max(nmap,maxval(pol_in(2)%vmap))
     allocate(wmap(nmap))
     wmap = .false. ! Initialization
     do isp = 1, 2
        do ii = 1, pol_in(isp)%ncond
           wmap(pol_in(isp)%cmap(ii)) = .true.
        enddo
        do ii = 1, pol_in(isp)%nval
           wmap(pol_in(isp)%vmap(ii)) = .true.
        enddo
     enddo
     ! W Gao dbg
     if(peinf%master) then
       write(*,'(a,i5)') " nmap = ", nmap
       write(*,'("   i    wmap(i)   ")')
       do ii = 1, nmap
          if(wmap(ii)) then
            write(*,'(i8," True")')  ii
          else
            write(*,'(i8," False")') ii
          endif
       enddo
     endif
  else
     allocate(wmap(1))
  endif
  if (min(pol_in(1)%ncond,pol_in(1)%nval,pol_in(2)%ncond,pol_in(2)%nval) < 0) then
     deallocate(wmap)
     allocate(wmap(1))
     nmap = 0
  endif

  !-------------------------------------------------------------------
  ! Read wave-function file.
  !
  init_gr = .false.
  if ( dft_code == PARATEC ) then
     call paratec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)
  else
     call parsec_wfn(gvec,kpt,nmap,nspin,wmap,init_gr)
  endif
  deallocate(wmap)
  if (trip_flag .and. nspin == 2) then
     write(lastwords,*) 'ERROR: cannot do TDLDA triplet kernel with ', &
          'spin polarized wavefunctions. Stop.'
     call die(lastwords)
  endif

  if(doisdf) then
     ! --- prepare some inputs for the ISDF method ---
     ! W Gao find the index of highest occupied orbital
     ihomo = 1
     do isp = 1, nspin
        do ikp = 1, kpt%nk
           do ii = 1, kpt%wfn(isp,ikp)%nstate
              if ( kpt%wfn(isp,ikp)%occ0(ii) > tol_occ .and. &
                   ihomo < ii) then
                 ihomo = ii
              endif
           enddo ! ii loop
        enddo ! ikp loop
     enddo ! isp loop
     ! if n_intp can not be found in rgwbs.in or invalid (i.e., less than the
     ! number of occupied states), then set it to a default value
     if(n_intp .lt. ihomo) then 
        n_intp = int(2.0 * ihomo)
     endif
     allocate(intp(n_intp))
     intp(1:n_intp) = 0
     ! --- find interpolation points for ISDF method ---
     call stopwatch(peinf%master, "before call cvt")
     call timacc(51,1,tsec)
     if(intp_type .eq. 1) then
        if (peinf%master) then
           write(*,*) " intp_type == 1"
           call cvt(gvec, kpt%rho, nspin, n_intp, intp)
        endif
     elseif(intp_type .eq. 2) then
        if (peinf%master) write(*,*) " intp_type == 2"
        call cvt_wfn(gvec, kpt%wfn, nspin, kpt%nk, n_intp, intp)
     else
        write(*,*) 'Type',intp_type,'method for finding interpolation points is',&
           ' not implememted so far. The default method will be used.'
     endif
     call MPI_BARRIER(peinf%comm,info)
     call timacc(51,2,tsec)
     call stopwatch(peinf%master, "after call cvt")
     ! pick out the interpolation points in reduced zone
     n_intp_r = 0
     do ipt = 1, n_intp
        if (mod(intp(ipt),gvec%syms%ntrans) == 1) n_intp_r = n_intp_r+1
     enddo
     allocate(intp_r(n_intp_r))
     ii = 0
     do ipt = 1, n_intp
        if (mod(intp(ipt),gvec%syms%ntrans) == 1) then
           ii = ii + 1
           intp_r(ii) = intp(ipt)
        endif
     enddo
     ! 
     if(peinf%master) write(*,*) " Finding interpolation points successfully. "
     ! broadcast intp to all processors 
     call MPI_bcast(intp(1),n_intp,MPI_INTEGER, peinf%masterid,peinf%comm,info)
     !
     ! ISDF will deal with all the pair products of wave functions as defined in
     ! pol_in(isp)%vmap(:) and pol_in(isp)%cmap(:)
     allocate(ncv( nspin,kpt%nk,gvec%syms%ntrans ))
     allocate(nv ( nspin,kpt%nk,gvec%syms%ntrans ))
     allocate(nc ( nspin,kpt%nk,gvec%syms%ntrans ))
     allocate(icv( gvec%syms%ntrans ))
     ! Initializing with zeros
     ncv = 0 
     nc  = 0 
     nv  = 0
     maxivv = 0 ! the absolute index of the highest occupied state
     maxicc = 0 ! the absolute index of the highest unoccupied state
     do isp = 1, nspin
        maxivv = max( maxval(pol_in(isp)%vmap(:)), maxivv)
        maxicc = max( maxval(pol_in(isp)%cmap(:)), maxicc)
     enddo
     allocate(pairmap(maxivv, maxicc, nspin, kpt%nk, gvec%syms%ntrans))
     do ikp = 1, kpt%nk
        do isp = 1, nspin
           do iv = 1, pol_in(isp)%nval
              ivv = pol_in(isp)%vmap(iv)
              irp = kpt%wfn(isp,ikp)%irep(ivv)
              nv(isp,ikp,irp) = nv(isp,ikp,irp) + 1
           enddo
           do ic = 1, pol_in(isp)%ncond
              icc = pol_in(isp)%cmap(ic)
              irp = kpt%wfn(isp,ikp)%irep(icc)
              nc(isp,ikp,irp) = nc(isp,ikp,irp) + 1
           enddo
           icv = 0
           do iv = 1, pol_in(isp)%nval
              do ic = 1, pol_in(isp)%ncond
                 ivv = pol_in(isp)%vmap(iv)
                 icc = pol_in(isp)%cmap(ic)
                 virp = kpt%wfn(isp,ikp)%irep(ivv)
                 cirp = kpt%wfn(isp,ikp)%irep(icc)
                 vcirp = gvec%syms%prod(virp, cirp)
                 icv(vcirp) = icv(vcirp) + 1
                 ! pairmap maps the real valence and conduction band index to 
                 ! the |vc> pair index obtained in isdf.f90 
                 pairmap(ivv,icc,isp,ikp,vcirp) = icv(vcirp)
                 ! invpairmap maps the 
           !      invpairmap(1,icv,isp,ikp) = ivv
           !      invpairmap(2,icv,isp,ikp) = icc
              enddo
           enddo
           ncv(isp, ikp, 1:gvec%syms%ntrans) = icv(1:gvec%syms%ntrans)
        enddo ! isp
     enddo ! ikp
     maxnc  = maxval( nc )
     maxnv  = maxval( nv )
     maxncv = maxval( ncv)
     write (6,*) " max Ncv", maxncv, " max Nv", maxnv, " max Nc", maxnc
     allocate(invpairmap(2, maxncv, nspin, kpt%nk, gvec%syms%ntrans))
     allocate(ivlist(maxnv, nspin, kpt%nk, gvec%syms%ntrans))
     allocate(iclist(maxnc, nspin, kpt%nk, gvec%syms%ntrans))
     ! Initialization with zeros
     pairmap    = 0     
     invpairmap = 0
     ivlist     = 0 
     iclist     = 0
     do ikp = 1, kpt%nk
        do isp = 1, nspin
           do iv = 1, pol_in(isp)%nval
              ivv = pol_in(isp)%vmap(iv)
              irp = kpt%wfn(isp,ikp)%irep(ivv)
              ivlist(iv,isp,ikp,irp) = ivv
           enddo
           do ic = 1, pol_in(isp)%ncond
              icc = pol_in(isp)%vmap(ic)
              irp = kpt%wfn(isp,ikp)%irep(icc)
              iclist(ic,isp,ikp,irp) = icc
           enddo
           icv = 0
           do iv = 1, pol_in(isp)%nval
              do ic = 1, pol_in(isp)%ncond
                 ivv = pol_in(isp)%vmap(iv)
                 icc = pol_in(isp)%cmap(ic)
                 virp = kpt%wfn(isp,ikp)%irep(ivv)
                 cirp = kpt%wfn(isp,ikp)%irep(icc)
                 vcirp = gvec%syms%prod(virp, cirp)
                 icv(vcirp) = icv(vcirp) + 1
                 invpairmap(1,icv(vcirp),isp,ikp,vcirp) = ivv
                 invpairmap(2,icv(vcirp),isp,ikp,vcirp) = icc
              enddo ! ic loop
           enddo ! iv loop
        enddo ! isp loop
     enddo ! ikp loop
  endif ! if (doisdf)

  !-------------------------------------------------------------------
  ! Calculate characters of representations.
  !
  if (kpt%lcplx) then
     call zcharac_group(gvec%syms,gvec,kpt,70,nspin,kpt%wfn(1,1)%nstate)
  else
     call dcharac_group(gvec%syms,gvec,kpt,70,nspin,kpt%wfn(1,1)%nstate)
  endif

  !-------------------------------------------------------------------
  ! Initialize pol and k_p structures.
  !
  allocate(pol(gvec%syms%ntrans,qpt%nk))
  pol(:,:)%ntr = 0
  allocate(k_p(gvec%syms%ntrans,qpt%nk))
  k_p(:,:)%ncol = 0
  k_p(:,:)%nbuff = nbuff
  k_p(:,:)%lcache = lcache
  k_p(:,:)%isdf =  doisdf
  ! Skip construction of kernel if we are doing RPA spectrum only.
  if (rpaonly) k_p(:,:)%ncol = -1

  call stopwatch(peinf%master,'Calling setup_g')
  call timacc(2,1,tsec)
  if (kpt%lcplx) then
     call zsetup_g(gvec,kpt,qpt,pol_in,pol,k_p,nspin,tdldacut,.false.)
  else
     call dsetup_g(gvec,kpt,qpt,pol_in,pol,k_p,nspin,tdldacut,.false.)
  call timacc(2,2,tsec)
  endif

  ! --- perform ISDF method to interpolate pair products of wave functions ---
  if(doisdf) then
     call stopwatch(peinf%master, "before call isdf")
     if (peinf%master) write(*,*) 'call isdf subroutine'
     kflag = 1 
     if ( trip_flag ) kflag = 3 
     if ( noxchange ) kflag = 2
     if ( nolda )     kflag = 0
     write(*, *) n_intp_r, maxncv, nspin, kpt%nk, gvec%syms%ntrans
     allocate(Cmtrx(n_intp_r, maxncv, nspin, kpt%nk, gvec%syms%ntrans)) 
     ! Cmtrx is intialized to zero in isdf_parallel.f90
     allocate(Mmtrx(n_intp_r, n_intp_r, nspin, nspin, kpt%nk, 2, gvec%syms%ntrans)) 
     ! Mmtrx is intialized to zero in isdf_parallel.f90
     verbose = .FALSE.
     call timacc(52,1,tsec)
     call isdf_parallel_sym(gvec, pol_in, kpt, n_intp_r, intp_r, &
       nspin, ncv, maxncv, &
       invpairmap, nv, maxnv, ivlist, nc, maxnc, iclist, &
       kflag, Cmtrx, Mmtrx, verbose)
     if(peinf%master) write(*,*) 'done isdf'
     call timacc(52,2,tsec)
     call stopwatch(peinf%master, "after call isdf")
  endif

  ! The outputs are Cmtrx and Mmtrx, which are used by k_integrate_isdf() for
  ! calculation of K(v,c,v',c') later !!
  ! --- finished ISDF ---

  !-------------------------------------------------------------------
  ! Print out warnings, information etc.
  !
  if (kpt%lcplx) tamm_d = .true.
  if (peinf%master) then
     write(6,'(/,a,/,/,a,/,2a,/)') repeat('-',65), &
          ' Polarizability input data: ', ' ', repeat('-',25)
     write(6,'(2a)') ' Number of transitions per representation ', &
          'in TDLDA polarizabillity:'
     write(6,'(8i8)') ((pol(irp,iq)%ntr,irp=1,gvec%syms%ntrans),iq=1,qpt%nk)
     write(6,'(a,i10,/)') ' total = ', sum(pol(:,:)%ntr)
     if (tdldacut > zero) then
        write(6,*) 'Energy cutoff applied in TDLDA polarizability = ', &
             tdldacut*ryd, ' eV'
     else
        write(6,*) 'No energy cutoff in TDLDA polarizability'
     endif
     if (nolda) write(6,'(/,a,/)') &
          ' LDA kernel is not included in polarizability'
     if (tamm_d) then
        write(6,'(2a,/)') ' Calculating TDLDA ', &
             'polarizability within the Tamm-Dancoff approximation.'
     else
        write(6,'(2a,/)') ' Not using the Tamm-Dancoff ', &
             'approximation in TDLDA polarizability.'
     endif
     !
     ! WARNING: for periodic boundary conditions, use a truncated Coulomb
     ! potential (see Hanke's article) in diagonalization. For that, set
     ! is_pbc = .false.
     ! In other applications (e.g., when calculating self-energies), one
     ! should not truncate the Coulomb potential.
     !
     write(6,'(a,/)') &
          ' Coulomb potential is being truncated in the TDLDA equation.'
     if (noxchange) write(6,'(a,/,a,/,a,/)') repeat('*',65), &
          ' Exchange kernel not included in TDLDA ',repeat('*',65)
     !
     ! Estimate memory usage.
     ! For diagonalization, we store 4 matrices: hamiltonian/eigenvectors, 
     ! temporary array (in eigensolver), Hartree kernel, and XC kernel.
     ! Parallelized diagonalization also uses a second temporary matrix.
     !
     mem1 = sum(kpt%wfn(:,:)%nmem)*two/real(nspin,dp) * &
          real(nspin*gvec%nr,dp)/two/131072.d0 / real(w_grp%npes,dp)
     if (kpt%lcplx) mem1 = mem1 * two
     write(6,'(/,a,f10.2,a)') ' Memory needed to store wavefunctions : ', &
           mem1,' MB/proc.'
     mem1 = real(5*gvec%nr,dp)/131072.d0
     if (kpt%lcplx) mem1 = mem1 * two
     write(6,'(a,a,f10.2,a)') ' Memory needed to calculate kernel ', &
          'matrix elements : ',mem1,' MB/proc.'
     xmax = 0
     do iq = 1, qpt%nk
        do irp = 1, gvec%syms%ntrans
           xsum = real(pol(irp,iq)%ntr * pol(irp,iq)%ntr,dp)
           if (xmax < xsum) xmax = xsum
        enddo
     enddo
     mem1 = xmax/1024.d0*3.d0/128.d0/r_grp%npes
     if (r_grp%npes > 1) mem1 = mem1*4.d0/3.d0
     write(6,'(/,a,f10.2,a/)') ' Memory needed for diagonalization : ', &
           mem1,' MB/proc.'
     write(6,'(a,/)') repeat('-',65)
     call flush(6)
  endif

  if (noxchange) trunc_c = .false.

  ii = 0
  xmax = zero
  call stopwatch(peinf%master, "after call calculate_tdlda")
  if (kpt%lcplx) then
     call zcalculate_tdlda(gvec,kpt,qpt,k_p,pol,nspin,ii, &
          tamm_d,nolda,rpaonly,trip_flag,noxchange,trunc_c,xsum,xmax, &
          Cmtrx, Mmtrx, n_intp_r, maxicc, maxivv, maxncv, pairmap)
  else
     call dcalculate_tdlda(gvec,kpt,qpt,k_p,pol,nspin,ii, &
          tamm_d,nolda,rpaonly,trip_flag,noxchange,trunc_c,xsum,xmax, &
          Cmtrx, Mmtrx, n_intp_r, maxicc, maxivv, maxncv, pairmap)
  endif
  call stopwatch(peinf%master, "after call calculate_tdlda")

  close(outdbg)

  ! Deallocate arrays for ISDF method
  if(doisdf) then
    deallocate(pairmap)
    deallocate(invpairmap)
    deallocate(ivlist)
    deallocate(iclist)
    deallocate(nv)
    deallocate(nc)
    deallocate(Cmtrx)
    deallocate(Mmtrx)
    deallocate(ncv)
  endif

  !-------------------------------------------------------------------
  ! Time accounting.
  !
  ii = 3
  jj = 7
  allocate(routnam(ii+jj))
  allocate(timerlist(ii+jj))
  routnam(1)='SETUP_T:'      ; timerlist(1)=2
  routnam(2)='KERNEL:'       ; timerlist(2)=3
  routnam(3)='DIAG_POL:'     ; timerlist(3)=4

  routnam(4)='POISSON_FFT:'         ; timerlist(4)=11
  routnam(5)='EIGENSOLVER:'         ; timerlist(5)=12
  routnam(6)='INTEGRATION:'         ; timerlist(6)=13
  routnam(7)='Find intp pts:'       ; timerlist(7)=51
  routnam(8)='ISDF_PARALLEL:'       ; timerlist(8)=52
  routnam(9)='Calc intp vectors:'   ; timerlist(9)=53
  routnam(10)='Calc <zeta|K|zeta>:' ; timerlist(10)=54


  call timacc(1,2,tsec)
  call finalize(peinf%master,peinf%comm,ii,jj,routnam,timerlist)

end program tdlda
!===================================================================
