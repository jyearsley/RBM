Module Block_Ice_Snow
!
! Module with ice and snow parameters and variables
!
! Integer variables 
!
!
! Logical variables 
!
    logical :: ICE(1000), SNOW(1000), SUB_ZERO(1000)
!
! Real variables
!
    real, parameter                 :: ice_albedo = 0.25
    real, parameter                 :: ice_cndctvy = 2.17   ! W/m/deg K  
    real, parameter                 :: h2o_dnsty = 1000.0  ! Kg/m**3
    real, parameter                 :: ice_dnsty = 917.0   ! Kg/m**3
    real, parameter                 :: ice_min_thck = 0.01 ! meters
    real, parameter                 :: sno_dnsty = 330.0   ! Kg/m**3
    real, parameter                 :: sens_h2o_ice = 0.280e-03 !Saucier et al
    real, parameter                 :: ltnt_ice_air = 0.840e-03 !Saucier et al
    real, parameter, dimension(2)   :: sens_ice_air = (/0.69e-03,1.12e-03/)       
!
    real, dimension(:,:,:), allocatable :: ice_temp, ice_thick
!
    real, parameter, dimension(2)   :: snow_albedo = (/0.85,0.60/)
    real, parameter, dimension(2)   :: ice_extnct  = (/1.5,20.0/) ! 1/meters
    real, parameter, dimension(2)   :: snow_extnct = (/6.0,20.0/) ! 1/meters
!	
end module Block_Ice_Snow
