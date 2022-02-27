module Block_Energy
!
!   Energy budget variables
!
!   Incoming short wave radiation, Watts/meter**2
!
    real, dimension(:), allocatable     :: qns
!
!   Incoming atmospheric radiation, Watts/meter**2
!
    real, dimension(:), allocatable     :: qna
!
!   Air temperature at surface, deg. C
!
    real, dimension(:), allocatable     :: dbt
!  
!   Wind speed, m/sec
!
    real, dimension(:), allocatable     :: wind
!
!   Vapor pressure of air at surface, kPa
!
    real, dimension(:), allocatable     :: ea
!
!   Air pressure at surface, kPa
!
    real, dimension(:), allocatable     :: press 
!
! Simulated water temperatures 
!   
    real, dimension(:),   allocatable   :: T_head,T_smth,T_trib 
    real, dimension(:,:), allocatable   :: Q_nps
    real, dimension(:,:), allocatable   :: temp_trib
    real, dimension(:,:), allocatable   :: temp_nps,thermal
    !
    real, dimension(:,:,:), allocatable :: temp

!
    real, dimension (:), allocatable   :: mu,alphamu,beta,gmma,smooth_param
!
    real, dimension (:,:), allocatable :: f1_tilde
!   Some important constants
!
      real             :: T_0
      real             :: kcal_Wsec = 4184.0
      real             :: lvp                    ! Latent heat of vaorization - J(W sec)/kg
      real,parameter   :: lvs = 2.834e06         ! Latent heat of sublimation                    
      real,parameter   :: rho = 1000.0           ! Kg/meter**3
      real,parameter   :: wind_fctr = 1.0        ! dimensionless
      real,parameter   :: rho_Cp = 1000.0*4184.0 ! (Kg/m**3)*(Joules/(Kg*degK)
      real,dimension(2),parameter :: evrate = (/1.5e-8,0.75e-8/) ! 1/kPa
!
end module Block_Energy  