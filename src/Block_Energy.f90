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
    real, dimension(:),allocatable :: press 
!
    real, dimension(:),allocatable :: mu,alphamu,beta,gmma,smooth_param

!   Some important constants
!
      real             :: lvp,rb,rho
      real             :: dlta,dlta1,dlta2,t0,t1,t2,x1,x2
      real,parameter   :: pi=3.14159
      real,parameter   :: rho_Cp = 1000.*4184.0 !rho*Cp(kg/meter**3)*(J/KG*degK)       
!
end module Block_Energy  
