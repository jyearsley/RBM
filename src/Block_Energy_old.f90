module Block_Energy
!
!   Energy budget variables
!
!   Incoming short wave radiation, kcal/m**2/sec
!
    real, dimension(:), allocatable::q_ns
!
!   Incoming atmospheric radiation, kcal/m**2/sec
!
    real, dimension(:), allocatable::q_na
!
!   Air temperature at surface, deg. C
!
    real, dimension(:), allocatable::dbt
!  
!   Wind speed, m/sec
!
    real, dimension(:), allocatable::wind
!
!   Vapor pressure of air at surface, mb
!
    real, dimension(:), allocatable::ea
!
!   Air pressure at surface, mb
!
    real, dimension(:), allocatable::press 

!
    real, dimension (:), allocatable::mu,alphamu,beta,gmma,smooth_param

!   Some important constants
!
      real             :: lvp                    ! Wikipedia - joules(W sec)/kg
      real             :: rb                 
      real             :: rho = 1000.0           ! Kg/meter**3
      real             :: wind_fctr = 1.0        ! dimensionless
      real,parameter   :: rho_Cp = 1000.0*4184.0 ! (Kg/m**3)*(Joules/(Kg*degK)
      real,dimension(2),parameter :: evrate = (/1.5e-8,0.75e-8/) ! 1/kPa
!
end module Block_Energy  