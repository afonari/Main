!-----------------------------------------------------------------------
subroutine dump_cart_disp (nat,ntyp,amass,ityp,dyn)
!-----------------------------------------------------------------------
  !
  !   diagonalise the dynamical matrix
  !   On input:  amass = masses, in amu
  !   On output: w2 = energies, z = displacements
  !
  use kinds, only: dp
  use constants, only: amu_ry, ry_to_cmm1, ry_to_thz
  implicit none
  ! input
  integer nat, ntyp, ityp(nat)
  complex(DP) dyn(3,3,nat,nat)
  real(DP) amass(ntyp)
  ! output
  ! NAN
  ! local
  real(DP) diff, dif1, difrel
  integer nat3, na, nta, ntb, nb, ipol, jpol, i, j
  complex(DP), allocatable :: dyn2(:,:), z(:,:)

  integer nu_i, nu_j, mu, iout
  real(DP) w1
  real(DP), allocatable :: u(:,:), m(:,:), w2(:)

  nat3 = 3*nat
  iout = 16
  open (unit=iout,file='dump_cart_disp',status='unknown',form='formatted')


  allocate(u(nat3, nat3))
  u(:,:) = 0.d0
  DO i = 1,nat3
    u(i,i) = 1.0d0
  ENDDO

  !  stolen from rigid.f90
  !  fill the two-indices dynamical matrix
  !
  allocate(dyn2 (nat3, nat3))
  !
  do na = 1,nat
     do nb = 1,nat
        do ipol = 1,3
           do jpol = 1,3
              dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = dyn(ipol,jpol,na,nb)
           end do
        end do
     end do
  end do
  !
  !  impose hermiticity
  !
  diff = 0.d0
  difrel=0.d0
  do i = 1,nat3
     dyn2(i,i) = CMPLX( DBLE(dyn2(i,i)),0.d0,kind=DP)
     do j = 1,i - 1
        dif1 = abs(dyn2(i,j)-CONJG(dyn2(j,i)))
        if ( dif1 > diff .and. &
             max ( abs(dyn2(i,j)), abs(dyn2(j,i))) > 1.0d-6) then
           diff = dif1
           difrel=diff / min ( abs(dyn2(i,j)), abs(dyn2(j,i)))
        end if
        dyn2(i,j) = 0.5d0* (dyn2(i,j)+CONJG(dyn2(j,i)))
        dyn2(j,i) = CONJG(dyn2(i,j))
     end do
  end do
  if ( diff > 1.d-6 ) write (iout,'(5x,"Max |d(i,j)-d*(j,i)| = ",f9.6,/,5x, &
       & "Max |d(i,j)-d*(j,i)|/|d(i,j)|: ",f8.4,"%")') diff, difrel*100
  !
  ! stolen from dyndiar.f90
  !
  !  fill the mass matrix (masses are in amu, amu_ry converts to a.u.)
  !
  allocate(m(nat3, nat3))
  DO nu_i = 1,nat3
     DO nu_j = 1,nat3
        m(nu_i,nu_j) = 0.0d0
        DO mu = 1,nat3
           na = (mu-1)/3+1
           m(nu_i,nu_j) = m(nu_i,nu_j) + amu_ry*amass(ityp(na))*u(mu,nu_i)*u(mu,nu_j)
        ENDDO
     ENDDO
  ENDDO

  allocate( z(nat3,nat3), w2(nat3) )
  CALL cdiaghg (nat3, nat3, dyn, m, nat3, w2, z)
  !  write frequencies
  !
  WRITE( iout,'(5x,"diagonalizing the dynamical matrix ..."//)')
  WRITE( iout,'(1x,74("*"))')
  !
!  dynout (:,:) = 0.0d0
!  DO nu_i = 1,nat3
!!     w1 = sqrt(abs(w2(nu_i)))
!     IF (w2(nu_i)<0.0) w1 = -w1
!     WRITE( stdout,9010) nu_i, w1*ry_to_thz, w1*ry_to_cmm1
     !  bring eigendisplacements in cartesian axis

     !SF
!     do na = 1,nat
!        write (stdout,9020) (z((na-1)*3+ipol,nu_i)/sqrt(amu_ry * amass(ityp(na))),ipol=1,3)
!     end do

     !do na = 1,nat
     !   write(stdout,*)
     !   write (stdout,9020) (u((na-1)*3+ipol,nu_i),ipol=1,3)
     ! end do
     ! end SF

!     DO mu = 1,3*nat
!        DO i = 1,nmodes
!           dynout(mu,nu_i) = dynout(mu,nu_i) + z(i,nu_i)*u(mu,i)
!        ENDDO
!     ENDDO
!  ENDDO
!  WRITE( stdout,'(1x,74("*"))')

  deallocate(m)
  deallocate(u)
  deallocate(dyn2)
  return
end subroutine dump_cart_disp



