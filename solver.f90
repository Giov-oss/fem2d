module solver
contains
  subroutine GAUSS(c,b,x,n)
    implicit none
    integer i,j,n
    real,dimension(n,n)::c
    real,dimension(n)::b
    real,dimension(n)::x
    real,dimension(n,n+1)::a
  a(:,1:n)=c
  a(:,n+1)=b
  do i=1,n
    a(i,:)=a(i,:)/a(i,i)
    do j=1,n
      if (i /= j) then
        a(j,:)=a(j,:)-a(i,:)*a(j,i)
      end if
    end do
  end do
  x(:)=a(:,n+1)
  end subroutine
end module
