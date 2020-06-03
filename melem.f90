module elemental
  contains
  subroutine matrizelemental(melem,kelem,fuente,e,nelementos,nnodos,nodos,elementos)
    use energia
    implicit none
    !matriz elemental
      integer :: e,nelementos,nnodos,i,l,m
      real :: nodos(nnodos,3), elementos(nelementos,5)
      real,dimension(4,4) :: n,dnk,dne
      real :: dn(2,4,4),coord(4,2),jac(2,2,4),ijac(2,2,4),dnxy(2,4,4)
      real :: melem(4,4),kelem(4,4),fuente(4,1)
      real :: djac(4),masa,kc,kd,f
  
      call shape_functions(n,dnk,dne)
      call coordenadas(e,coord,nnodos,nodos,nelementos,elementos)

    !jacobiano e integración
      do i=1,4
        dn(1,:,i)=dnk(:,i)
        dn(2,:,i)=dne(:,i)
        jac(:,:,i)=matmul(dn(:,:,i),coord)
        djac(i)=jac(1,1,i)*jac(2,2,i)-jac(2,1,i)*jac(1,2,i)
        ijac(1,1,i)=jac(2,2,i)
        ijac(1,2,i)=-jac(1,2,i)
        ijac(2,1,i)=-jac(2,1,i)
        ijac(2,2,i)=jac(1,1,i)
        ijac(:,:,i)=ijac(:,:,i)/djac(i)
        dnxy(:,:,i)=matmul(ijac(:,:,i),dn(:,:,i))
      end do
  
    do l=1,4
      f=0.
      do m=1,4
        masa=0.
        kd=0.
        kc=0.
        do i=1,4
          masa=masa+n(l,i)*n(m,i)*djac(i)
          kd=kd+(k/(rho*c))*(dnxy(1,l,i)*dnxy(1,m,i)+dnxy(2,l,i)*dnxy(2,m,i))*djac(i)
          kc=kc+n(l,i)*(u*dnxy(1,m,i)+v*dnxy(2,m,i))*djac(i)
        end do
        kelem(l,m)=difus*kd+convec*kc
        melem(l,m)=time*masa
      end do
      do i=1,4
        f=f+(q/(rho*c))*n(l,i)*djac(i)
      end do
      fuente(l,1)=f
    end do
  
  end subroutine matrizelemental

  subroutine shape_functions(n,dnk,dne)
    implicit none 
    integer :: i,j
    integer,parameter :: ngp=4
    real,dimension(ngp,2) :: gqke=reshape((/-sqrt(1./3.),sqrt(1./3.),-sqrt(1./3.),&
    sqrt(1./3.),-sqrt(1./3.),-sqrt(1./3.),sqrt(1./3.),sqrt(1./3.)/),(/4,2/))
    integer::ksii,etai
    real,dimension(4,4) :: n,dnk,dne
    do j=1,ngp
      do i=1,4
        select case (i)
          case (1)
          ksii=-1
          etai=-1
          case (2)
          ksii=1
          etai=-1
          case (3)
          ksii=1
          etai=1
          case (4)
          ksii=-1
          etai=1
        end select
      n(i,j)=0.25*(1.+ksii*gqke(j,1))*(1.+etai*gqke(j,2))
      dnk(i,j)=0.25*ksii*(1+etai*gqke(j,2))
      dne(i,j)=0.25*etai*(1+ksii*gqke(j,1))
      end do
    end do
  end subroutine
  
  subroutine coordenadas(e,coord,nnodos,nodos,nelementos,elementos)
    implicit none
    integer :: j,e,n1,n2,n3,n4,nnodos,nelementos
    real :: coord(4,2),nodos(nnodos,3),elementos(nelementos,5)
    n1=int(elementos(e,2))
    n2=int(elementos(e,3))
    n3=int(elementos(e,5))
    n4=int(elementos(e,4))
    do j=1,nnodos
      if (int(nodos(j,1))==n1)      then
        coord(1,:)=nodos(j,2:3)
      else if (int(nodos(j,1))==n2) then
        coord(2,:)=nodos(j,2:3)
      else if (int(nodos(j,1))==n3) then
        coord(3,:)=nodos(j,2:3)
      else if (int(nodos(j,1))==n4) then
        coord(4,:)=nodos(j,2:3)
      end if
    end do
  end subroutine
  
  subroutine nlocaltoglobal(elementos,nelementos,i,eg)
    implicit none
    integer :: i,n1,n2,n3,n4,nelementos
    real,dimension(nelementos,5) :: elementos
    integer,dimension(4) :: eg
      n1=int(elementos(i,2))
      n2=int(elementos(i,3))
      n3=int(elementos(i,5))
      n4=int(elementos(i,4))
      eg=(/n1,n2,n3,n4/)
  end subroutine
  
end module
  
  
  