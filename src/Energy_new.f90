      SUBROUTINE ENERGY
     &           (Tsurf,Qsurf,wind_fctr,A,B,ncell)
!
    use Block_Energy
   implicit none
   integer::i,ncell,nd
   real                :: A,B,e0,Qsurf,QCONV,QEVAP,QWS
   real                :: td,Tsurf
   real, dimension(2)  :: q_fit, T_fit
!     
      td=nd
!
! Testing effect of hysteresis with output as Farmington_3  4/2/2018
!
      evap_rate=evrate(1)
      if (nd > 180) evap_rate=evrate(2)
      T_fit(1) = Tsurf-1.0
      T_fit(2) = Tsurf+1.0
      if (T_fit(1) .lt. 0.50) T_fit(1) = 0.50
      if (T_fit(2) .lt. 0.50) T_fit(2) = 1.00
      do i=1,2
         T_kelvin = T_fit(i) + 273.0
!
! Vapor pressure at water surface - kPa (Tetens)
         e0=0.61078*exp((17.27*T_fit(i))/T_Kelvin)
!
! Bowen ratio
         rb = 0.491*((dbt(nd)-T_fit(i))*press(ncell)/
!
! Latent heat of vaporization - joules (Wsec)/kg
         lvp = kcal_Wsec * (597.0 - T_fit(i))
!
! Evaporative heat flux
         QEVAP=wind_fctr*rho*lvp*evap_rate*WIND(ncell)
         if(QEVAP.lt.0.0) QEWVAP=0.0
!
! Convective heat flux
         QCONV=rb*QEVAP
         QEVAP=QEVAP*(E0-EA(ncell))
! 
! Back radiation from the water surface - W/m**2
         QWS=280.23+6.1589*T_fit(i)
!
! Thermal energy budget for i = 1,2
         q_fit(i)=QNS(ncell)+0.97*QNA(ncell)-QWS-QEVAP+QCONV
      end do
c
! Linear relationship for equilibrium temperature and rate constant
!
      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))
     .     /(T_fit(1)-T_fit(2))
!
      Qsurf=0.5*(q_fit(1)+q_fit(2))
      RETURN
      END