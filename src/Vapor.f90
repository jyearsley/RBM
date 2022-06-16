REAL FUNCTION VAPOR
use Block_Ice_Snow
use Block_Energy
!
Implicit None
!
q10 = eps_air_h2o*ea(ncell)/(press(ncell)-comp_eps*ea(ncell))
q_ice = eps_air_h2o*e0/(press(ncell)-comp_eps*e0) 
END REAL FUNCTION VAPOR
