!------------------------------------------------------------------------------
!------------------------------------------------------------------main program
program fem_galerkin
  use energia
  use domain
  use elemental
  use solver
  implicit none
  !discretización
  real,allocatable,dimension(:,:):: nodos,elementos,bc
  !matrices
  real,allocatable,dimension(:,:) :: mglobal,dmglobal,kglobal,dkglobal,rhs
  integer :: alloc_stat
  integer :: nnodos,nelementos,nbc=0,j,i,h,ts=(tf-t0)/dt
  real :: melem(4,4),kelem(4,4),fuente(4,1) !matriz elemental
  integer,dimension(4) :: eg !nodos del elemento i
  !matrices sist. ecs.
  real,allocatable,dimension(:,:) :: temp,ml,mr,mlr,mrr !matriz reducida
  real,allocatable,dimension(:) :: rhsr,var
  integer,allocatable,dimension(:) :: reductor !operador reductor
  real :: coord(4,2)
!------------------------------------------------------------------------------
!--------------------------------------------------------------lectura de datos
  call grid()
  call getdata('data/nods.dat',nodos,3)
  call getdata('data/elmnts.dat',elementos,5)
  call getdata('data/bc.dat',bc,6)

  nnodos=size(nodos)/3
  nelementos=size(elementos)/5

  allocate (mglobal(nnodos, nnodos), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"
  mglobal=0.

  allocate (dmglobal(nnodos, nnodos), stat = alloc_stat)
  if (alloc_stat /= 0) stop "** * not enough memory ***"
  dmglobal=0.

  allocate (kglobal(nnodos, nnodos), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"
  kglobal=0.

  allocate (dkglobal(nnodos, nnodos), stat = alloc_stat)
  if (alloc_stat /= 0) stop "** * not enough memory ***"
  dkglobal=0.

  allocate (rhs(nnodos, 1), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"
  rhs=0.
!------------------------------------------------------------------------------
!--------------------------------------------------------------matriz elemental
  do i=1,nelementos
    call matrizelemental(melem,kelem,fuente,i,nelementos,nnodos,nodos,elementos)
!------------------------------------------------------------------------------
!-----------------------------------------------------------------matriz global
    call nlocaltoglobal(elementos,nelementos,i,eg)
    do j=1,4
      do h=1,4
        dmglobal(eg(j),eg(h))=melem(j,h)
        dkglobal(eg(j),eg(h))=kelem(j,h)
        mglobal=mglobal+dmglobal
        kglobal=kglobal+dkglobal
        dmglobal=0.
       dkglobal=0.
      end do
      rhs(eg(j),1)=rhs(eg(j),1)+fuente(j,1)
    end do
  end do
  !call imprimir_matriz(mglobal,nnodos,nnodos)
  !call imprimir_matriz(kglobal,nnodos,nnodos)
!-----------------------------------------------------------------------------
!------------------------------------------------cond. de contorno e iniciales
  if (time==0) then
    ts=1
    theta=1
  end if

  allocate (temp(nnodos,ts+1), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"

  temp(nnodos,1)=0.
  do i=1,nnodos
    if (bc(i,2)==1) then
      temp(i,:)=bc(i,3)
      nbc=nbc+1
    else if (bc(i,2)==2) then
      rhs(i,1)=rhs(i,1)+bc(i,4)
    else if (bc(i,2)==3) then
      rhs(i,1)=rhs(i,1)+bc(i,6)
      kglobal(i,i)=kglobal(i,i)-bc(i,5)
    end if
  end do

  allocate (reductor(nnodos-nbc), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"

  j=1
  do i=1,nnodos
    if (bc(i,2)/=1) then
      reductor(j)=int(bc(i,1))
      j=j+1
    end if
  end do
!-----------------------------------------------------------------------------
!--------------------------------------------------------------sistema de ecs.
  allocate (ml(nnodos,nnodos), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"
  allocate (mr(nnodos,nnodos), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"
  allocate (mlr(nnodos-nbc,nnodos-nbc), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"
  allocate (mrr(nnodos-nbc,nnodos-nbc), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"
  allocate (rhsr(nnodos-nbc), stat = alloc_stat)
  if (alloc_stat /= 0) stop "*** not enough memory ***"

  ml=mglobal/dt+theta*kglobal
  mr=-mglobal/dt+(1-theta)*kglobal
  mlr=ml(reductor,reductor)
  mrr=mr(reductor,reductor)

 allocate (var(nnodos-nbc), stat = alloc_stat)
 if (alloc_stat /= 0) stop "*** not enough memory ***"

  do i=2,ts+1
    rhsr=rhs(reductor,1)-matmul(mrr,temp(reductor,i-1))
    do j=1,nnodos
      if (bc(j,2)==1) then
       rhsr=rhsr-(ml(reductor,int(bc(j,1)))*bc(j,3))-(mr(reductor,int(bc(j,1)))*bc(j,3))
      end if
    end do
    !call imprimir_matriz(rhsr,nnodos-nbc,1)
    call gauss(mlr,rhsr,var,nnodos-nbc)
    !temp(int(bc(:,1)),i)=bc(:,2)
    temp(reductor,i)=var
  end do
 !-----------------------------------------------------------------------------
 !--------------------------------------------------------------------resultado
  open(20,file='data/res.dat',status='unknown',action='write')
   do i=1,nnodos
     write(20,*)i,nodos(i,2),nodos(i,3),temp(i,2:)
   end do
  close(20)
  open(21,file='data/resml.dat',status='unknown',action='write')
  do i=1,nelementos
    call coordenadas(i,coord,nnodos,nodos,nelementos,elementos)
    do j=1,4
      write(21,*)int(elementos(i,j+1)),coord(j,1),coord(j,2),temp(int(elementos(i,j+1)),ts+1)
    end do
  end do
  close(21)
 !-----------------------------------------------------------------------------
 !-------------------------------------------------------------fin del programa
 deallocate (nodos)
 deallocate (elementos)
 deallocate (mglobal)
 deallocate (dmglobal)
 deallocate (kglobal)
 deallocate (dkglobal)
 deallocate (ml)
 deallocate (mlr)
 deallocate (mr)
 deallocate (mrr)
 deallocate (rhs)
 deallocate (rhsr)
 deallocate (var)
 deallocate (temp)
 deallocate (reductor)
 deallocate (bc)
!------------------------------------------------------------------------------
end program
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!-------------------------------------------------------------------subroutines
subroutine imprimir_matriz(matriz,fil,col)
implicit none
integer i,fil,col
real matriz(fil,col)
write(*,*) "----------------------------------------------------"
do i=1,fil
  write (*,'(101f9.4)') matriz(i,:)
end do
write(*,*) "----------------------------------------------------"
end subroutine