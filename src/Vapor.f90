SUBROUTINE VAPOR(q_ice,q10,T_p,ncell,ns)
use Block_Ice_Snow
use Block_Energy
!
Implicit None

integer             :: ncell,ns
!
real                :: ea_pascal,e0,expnt,pa_pascal,q10,q_ice,T_p
!
ea_pascal = 1000.*ea(ncell)
pa_pascal = 1000.*press(ncell)
q10 = eps_air_h2o*ea_pascal/(pa_pascal-comp_eps*ea_pascal)
!
expnt = a_vapor(1)*(T_p-273.160)/(T_p-b_vapor(1))
e0 = 611.*(10.**expnt)
q_ice = eps_air_h2o*e0/(pa_pascal-comp_eps*e0) 
!
END SUBROUTINE VAPOR
