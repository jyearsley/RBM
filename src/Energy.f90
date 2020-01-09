      SUBROUTINE ENERGY (Tsurf,qsurf,wind_fctr,A,B,ncell)
!
      use Block_Energy
      implicit none
      integer                       :: i,ncell,nd
      real                          :: Tsurf,Qsurf,wind_fctr,A,B
      real                          :: T_kelvin
      real, dimension(2)            :: q_fit,T_fit
      real                          :: evap_rate
      real                          :: e0,q_evap,q_conv,q_ws
      real, dimension(2), parameter :: evrte =(/1.5e-11,0.75e-11/)
      real, parameter               :: RHO=1000.,EVRATE=1.5e-11
!
! Testing effect of hysteresis
!
      evap_rate=evrte(1)
      if (nd > 180) evap_rate=evrte(2)
      T_fit(1) = Tsurf-1.0
      T_fit(2) = Tsurf+1.0
      if (T_fit(1) .lt. 0.50) T_fit(1) = 0.50
      if (T_fit(2) .lt. 0.50) T_fit(2) = 1.00
      do i=1,2
         T_kelvin = T_fit(i) + 273.0
!
! Vapor pressure at water surface
         e0=2.1718E10*EXP(-4157.0/(T_kelvin-33.91))
!
! Bowen ratio
         rb=PF*(DBT(ncell)-T_fit(i))
!
! Latent heat of vaporization
         lvp=1.91846e06*(T_kelvin/(T_kelvin-33.91))**2
!
! Evaporative heat flux
         Q_EVAP=wind_fctr*rho*lvp*evap_rate*WIND(ncell)
         if(q_evap.lt.0.0) q_evap=0.0
!
! Convective heat flux
         Q_CONV=rb*QEVAP
         Q_EVAP=Q_EVAP*(E0-EA(ncell))
! 
! Back radiation from the water surface
         Q_WS=280.23+6.1589*T_fit(i)
!
! Thermal energy budget for i = 1,2
         q_fit(i)=Q_NS(ncell)+0.97*Q_NA(ncell)-Q_WS-Q_EVAP+Q_CONV
      end do
!
! Linear relationship for equilibrium temperature and rate constant
!
      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))                  &
       /(T_fit(1)-T_fit(2))
!
      qurf=0.5*(q_fit(1)+q_fit(2))
      RETURN
      END