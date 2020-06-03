!------------------------------------------------------------------------------
!-----------------------------------------------------------ecuación de energía
module energia
  implicit none
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------master
  integer,parameter :: time=0
  integer,parameter :: convec=1
  integer,parameter :: difus=1
  integer,parameter :: supg=0
  !----------------------------------------------------------------------------
  !--------------------------------------------------------------------temporal
  real,parameter :: t0=0.
  real,parameter :: tf=3.
  real,parameter :: dt=0.1
  real :: theta=0.5
  !----------------------------------------------------------------------------
  !-----------------------------------------------------------------propiedades
  real,parameter :: rho=1.
  real,parameter :: c=1.
  real,parameter :: k=1.
  real,parameter :: u=10.,v=0.
  real,parameter :: q=0.
  !----------------------------------------------------------------------------
  !--------------------------------------------------------------discretización
  integer,parameter :: xelem=10
  integer,parameter :: yelem=10
  real,parameter :: lx=1.
  real,parameter :: ly=1.

end module energia