module stale
  integer(4) :: n, &      ! cluster size n x n
                nheavy, & ! number of heavy particles  
                n_spec    ! number of species
  real(8) :: temperature, U, cp
  real(8), parameter :: pi=3.14159265359_8
end module stale
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
module matrices
  use stale
  complex(8),dimension(:,:),allocatable :: ham
  integer(1),dimension(:,:),allocatable :: hconf
  integer(4) :: ix_from,ix_to,iy_from,iy_to ! hoppings
end module matrices

program fk !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use stale
use matrices
use, intrinsic :: ISO_FORTRAN_ENV
implicit none
!
integer(4) :: i,nn,j,ix,iy,k,iomega,mm,itherm
real(8),dimension(:),allocatable :: ens,rho,rho0
integer(8), allocatable :: m(:)
real(8) :: free_en,old_free_en,free_energy,x,l_number,energy,an,eps,omega,omax
character :: c
external energy
!
namelist /input/ n_spec,n,nheavy,U,temperature,cp
open(10,file='fk.inp',status='OLD')
read(10,nml=input)
close(10)
!
cp=U/2. ! HALF-FILLING
nn=n*n
mm=1000 ! 2*mm+1 = number of omegas
allocate(ham(nn,nn),hconf(nn,nn),ens(nn),m(n),rho(-mm:mm),rho0(-mm:mm))
!
call system("echo $PPID > PID_RHO")
call random_seed

rho = 0._8
eps = 5.e-3
omax = 5._8
itherm = 100

hconf = 0
do ix = 1,n
  do iy = 1,n
    if (mod(ix+iy,2) == 0) hconf(ix,iy) = 1
  enddo
enddo
rho0 = 0.0_8 ! DOS for the perfect checkerboard
do k=1,500
  call init_hamiltonian(.true.)
  call ham_diag(ens)
  do iomega = -mm,mm
      omega = omax*iomega/mm
      do i = 1,nn
          rho0(iomega) = rho0(iomega) + 1._8/((ens(i) - cp - omega)**2 + eps**2)
      enddo
  enddo
enddo
rho0 = rho0 * eps/(2*pi) / (500*nn)

!
j = 0
open(10,file="configs.dat")
!!!!!!!!!!$OMP PARALLEL DO DEFAULT(PRIVATE)
mainloop: do
    do ix = 1,n
        read(10,*,end=99) m(ix)
    enddo
    read(10,*) c
    j = j + 1
    if (j > itherm) then
        hconf = 0
        do ix = 1,n
            do iy = 1,n
                k = n*(ix-1)+iy
                if (btest(m(ix), iy-1)) then    
!                    ham(k,k) = U
                  hconf(ix,iy) = 1
                else
!                    ham(k,k) = 0.0_8
                  hconf(ix,iy) = 0
                endif                   
            enddo
        enddo
        call init_hamiltonian(.true.)
        call ham_diag(ens)
        do iomega = -mm,mm
            omega = omax*iomega/mm
            do i = 1,nn
                rho(iomega) = rho(iomega) + 1._8/((ens(i) - cp - omega)**2 + eps**2)
            enddo
        enddo
    endif
enddo mainloop
!!!!!!!!!!!!!$OMP END PARALLEL DO
99 continue
close(10)
rho = rho * eps/(2*pi)
rho = rho/((j - itherm)*nn)
open(11,file='rho.dat')
do iomega = -mm,mm
    write(11,*) omax*iomega/mm,rho(iomega),rho0(iomega)
enddo
close(11)
end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
subroutine init_conf
use stale
use matrices
implicit none
integer(4) :: i,ix,iy
real(8) :: x
open (14,file='tab_dom.dat',status='old',err=99)
do ix=1,n
  read (14,'(<n>I1)') (hconf(ix,iy),iy=1,n)
enddo
close(14)
return
!
99 hconf=0
do i=1,nheavy
  do
     call random_number(x); ix=int(n*x)+1
     call random_number(x); iy=int(n*x)+1
     if (hconf(ix,iy)==0) then
       hconf(ix,iy)=1
       exit
     endif
  enddo
enddo
end subroutine init_conf
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
subroutine init_hamiltonian(RBC)
use stale
use matrices
use, intrinsic :: ISO_FORTRAN_ENV
implicit none
!
integer(4) :: ix,iy,k1,k2
real(8) :: t=-1., & ! hopping integral
           phi
complex(8) :: cx,cy
logical :: RBC ! Random Boundary Conditions

if (RBC) then
  call random_number(phi)
  phi = 2*pi*phi
  cx = cmplx(cos(phi),sin(phi),kind=8)
  call random_number(phi)
  phi = 2*pi*phi
  cy = cmplx(cos(phi),sin(phi),kind=8)
else
  cx = cmplx(1._8,0._8,kind=real64)
  cy = cmplx(1._8,0._8,kind=real64)
endif

ham=cmplx(0._8,0._8,kind=real64)
do iy=1,n
  do ix=1,n
     k1=n*(ix-1)+iy
     if (ix < n) then
        k2=n*(ix+1-1)+iy
        ham(k1,k2)=t
        ham(k2,k1)=t
      else
        k2=iy
        ham(k1,k2)=t*cx
        ham(k2,k1)=t*conjg(cx)
      endif
     if (iy > 1) then
        k2=n*(ix-1)+iy-1
        ham(k1,k2)=t
        ham(k2,k1)=t   
     else
        k2=n*(ix-1)+n
        ham(k1,k2)=t*cy
        ham(k2,k1)=t*conjg(cy)
     endif
    if (hconf(ix,iy)==1) ham(k1,k1)=U
  enddo
enddo
end subroutine init_hamiltonian
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
subroutine hop
use stale
use matrices
real(8) :: x
do
   call random_number(x); ix_from=int(n*x+1)
   call random_number(x); iy_from=int(n*x+1)
   if (hconf(ix_from,iy_from)==1) exit
enddo
hconf(ix_from,iy_from)=hconf(ix_from,iy_from)-1
ham(n*(ix_from-1)+iy_from, n*(ix_from-1)+iy_from)= &
ham(n*(ix_from-1)+iy_from, n*(ix_from-1)+iy_from)-U
do
   call random_number(x); ix_to=int(n*x+1)
   call random_number(x); iy_to=int(n*x+1)
   if (hconf(ix_to,iy_to)==0) exit
enddo
hconf(ix_to,iy_to)=hconf(ix_to,iy_to)+1
ham(n*(ix_to-1)+iy_to, n*(ix_to-1)+iy_to)= &
ham(n*(ix_to-1)+iy_to, n*(ix_to-1)+iy_to)+U
end subroutine hop
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
subroutine ham_diag(ens) !Hamiltonian diagonalisation
use stale
use matrices
use blas95
use lapack95
implicit none
real(8),dimension(n*n) :: ens
complex(8),dimension(n*n,n*n) :: h
h=ham
!call syevd(h,ens)
call heevd(h,ens)
end subroutine ham_diag
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
subroutine save_results(free_en,en1,en2,an)
use stale
use matrices
implicit none
integer :: ix,iy
integer(8) :: m
real(8) :: x,free_en,en1,en2,an
open(14,file='tab_dom.dat')
do ix=1,n
  write(14,'(<n>I1)') (hconf(ix,iy),iy=1,n)
enddo
close(14)
open(15,file='free_en.dat',position='APPEND')
write(15,*) free_en/n_spec+cp*an,en1,en2
close(15)
open(15,file='n.dat',position='APPEND')
write(15,*) an
close(15)
open(16,file='configs.dat',position='APPEND')
do ix=1,n
   m=0
   do iy=1,n
      if (hconf(ix,iy)==1) m=ibset(m,iy-1)
   enddo
   write(16,*) m
enddo
write(16,*) '-'
close(16)
end subroutine save_results
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
subroutine revert_ham
use matrices
hconf(ix_from,iy_from)=hconf(ix_from,iy_from)+1
ham(n*(ix_from-1)+iy_from, n*(ix_from-1)+iy_from)= &
ham(n*(ix_from-1)+iy_from, n*(ix_from-1)+iy_from)+U
!
hconf(ix_to,iy_to)=hconf(ix_to,iy_to)-1
ham(n*(ix_to-1)+iy_to, n*(ix_to-1)+iy_to)= &
ham(n*(ix_to-1)+iy_to, n*(ix_to-1)+iy_to)-U
end subroutine revert_ham
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
function free_energy(ens)
use stale
implicit none
real(8),dimension(n*n) :: ens
real(8) :: free_energy
real(16) :: fe,expon
integer :: i
fe=0._16
!$OMP PARALLEL
!$OMP DO
do i=1,n*n
  expon=-(ens(i)-cp)/temperature
  if (abs(expon)<1.e-50_16) then
     fe=fe+log(2._16)
  else if (expon > 20) then
    fe=fe+expon
  else 
    fe=fe+log(1._16+exp(expon))
  endif
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
free_energy=-n_spec*temperature*fe
end function free_energy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
function energy(ens)
use stale
implicit none
real(8),dimension(n*n) :: ens
real(8) :: energy,fermi
integer :: i
energy=0._8
do i=1,n*n
  energy = energy + ens(i)*fermi(ens(i))
enddo
end function energy
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
function l_number(ens)
use stale
real(8),dimension(n*n) :: ens
real(8) :: l_number!,fermi
integer :: i
l_number=0._8
!$OMP PARALLEL
!$OMP DO
do i=1,n*n
  l_number=l_number+fermi(ens(i))
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
function fermi(e)
use stale
implicit none
real(8) :: fermi,e,ee
if (e==cp) then
  fermi=0.5
  return
endif
ee=(e-cp)/temperature
if (ee<-10) then
  fermi=1._8
else if (ee>10) then
  fermi=0._8
else 
  fermi=1._8/(exp(ee)+1)
endif
end function fermi
