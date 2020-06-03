module domain
  contains
  subroutine grid()
    use energia
    implicit none
    integer :: i,j
    integer,parameter :: nels=xelem*yelem
    integer,parameter :: nnodos=(xelem+1)*(yelem+1)
    real :: coord(nnodos,2),cx(xelem+1),cy(yelem+1),dx,dy
    real :: elements(nels,4)
    type edgebc
      integer :: bctype !(1: Dirichlet, 2: Neumann; 3: Robin)
      real :: pv,sv,beta,gamma !(pv: T, sv: q, beta,gamma: T+q)
    end type
    type(edgebc) :: cc(4)
    type(edgebc) :: corners(4)
    type bc
      integer :: nodo
      integer :: bctype
      real :: pv,sv,beta,gamma
    end type
    type(bc) :: bcs(nnodos)
    !-  
    cc(1)=edgebc(2,0,0,0,0) !pared inferior
    cc(2)=edgebc(1,1,0,0,0) !pared lateral derecha
    cc(3)=edgebc(2,0,0,0,0) !pared superior
    cc(4)=edgebc(1,0,0,0,0) !pared lateral izquierda
    !-
    corners(1)=cc(4) !esquina inferior izquierda
    corners(2)=cc(2) !esquina inferior derecha
    corners(3)=cc(2) !esquina superior derecha
    corners(4)=cc(4) !esquina superior izquierda
    !-
    !mallado: coordenadas, elementos y condiciones de borde  
    !-
    !coordenadas nodos
    dx=lx/xelem !tamaño de elemento en x
    dy=ly/yelem !tamaño de elemento en y
  
    cx=(/(dx*(i-1),i=1,xelem+1)/)
    cy=(/(dy*(i-1),i=1,yelem+1)/)

    do i=1,xelem+1
      do j=1,yelem+1
        coord(i+(xelem+1)*(j-1),1)=cx(i)
        coord(i+(xelem+1)*(j-1),2)=cy(j)
      end do
    end do
    !-
    !elementos
    do j=1,yelem
      do i=1,xelem
        elements(i+xelem*(j-1),1)=i+(xelem+1)*(j-1)
        elements(i+xelem*(j-1),2)=i+1+(xelem+1)*(j-1)
        elements(i+xelem*(j-1),3)=i+(xelem+1)*(j)
        elements(i+xelem*(j-1),4)=i+1+(xelem+1)*(j)
      end do
    end do
    !-
    !condiciones de borde: salida
    do i=1,nnodos
      if (coord(i,2)==0.) then
        bcs(i)=bc(i,cc(1)%bctype,cc(1)%pv,cc(1)%sv,cc(1)%beta,cc(1)%gamma)
      else if (coord(i,1)==lx) then
        bcs(i)=bc(i,cc(2)%bctype,cc(2)%pv,cc(2)%sv,cc(2)%beta,cc(2)%gamma)
      else if (coord(i,2)==ly) then
        bcs(i)=bc(i,cc(3)%bctype,cc(3)%pv,cc(3)%sv,cc(3)%beta,cc(3)%gamma)
      else if (coord(i,1)==0.) then
        bcs(i)=bc(i,cc(4)%bctype,cc(4)%pv,cc(4)%sv,cc(4)%beta,cc(4)%gamma)
      else
        bcs(i)=bc(i,0,0,0,0,0)
      end if
    end do
    !esquinas
    bcs(1)=bc(1,corners(1)%bctype,corners(1)%pv,corners(1)%sv,&
    corners(1)%beta,corners(1)%gamma)
    bcs(xelem+1)=bc(xelem+1,corners(2)%bctype,corners(2)%pv,corners(2)%sv,&
    corners(2)%beta,corners(2)%gamma)
    bcs(nnodos)=bc(nnodos,corners(3)%bctype,corners(3)%pv,corners(3)%sv,&
    corners(3)%beta,corners(3)%gamma)
    bcs(nnodos-yelem)=bc(nnodos-yelem,corners(4)%bctype,corners(4)%pv,corners(4)%sv,&
    corners(4)%beta,corners(4)%gamma)
    !-
    !-archivos
    open(69,file='data/elmnts.dat',status='unknown',action='write')
    open(70,file='data/nods.dat',status='unknown',action='write')
    open(71,file='data/bc.dat',status='unknown',action='write')
  
    do i=1,nels
      write(69,*)i,elements(i,:)
    end do

    do i=1,nnodos
      write(70,*)i,coord(i,:)
    end do
    do i=1,nnodos
      !if (bcs(i)%bctype/='nap') then
        write(71,*)bcs(i)
      !end if
    end do
  
    close(69)
    close(70)
    close(71)
    
    end subroutine grid

    subroutine getdata(name,matriz,c)
      implicit none
      integer,parameter::n=25921
      character(len=*):: name
      integer :: i,j,alloc_stat,c
      real,allocatable,dimension(:,:) :: dummy,matriz
      open(69,file=name,status='old',action='read')
      allocate (dummy(n, c), stat = alloc_stat)
      if (alloc_stat /= 0) stop "*** not enough memory ***"
      j=0
      do i=1,n
        if (c==3) then
          read(69,*,end=25)dummy(i,1),dummy(i,2),dummy(i,3)
        else if (c==5) then
          read(69,*,end=25)dummy(i,1),dummy(i,2),dummy(i,3),dummy(i,4),dummy(i,5)
        else if (c==6) then
          read(69,*,end=25)dummy(i,1),dummy(i,2),dummy(i,3),dummy(i,4),dummy(i,5),dummy(i,6)
        end if
      j=j+1
      end do
      allocate (matriz(j, c), stat = alloc_stat)
      if (alloc_stat /= 0) stop "*** not enough memory ***"
      25 matriz=dummy(1:j,1:c)
      close(69)
    end subroutine getdata

  end module
  
  
  