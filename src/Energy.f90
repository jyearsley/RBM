      SUBROUTINE ENERGY   (Tsurf, qsurf, nd, ncell, ICE)
!
    use Block_Energy
   implicit none
   integer             ::i, ncell, nd
   logical             :: Ice
   real                :: A, B, e0, evap_rate, qsurf, QCONV, QEVAP, QWS
   real                :: T_Kelvin, T_rb, Tsurf, vpr_diff
   real, dimension(2)  :: q_fit, T_fit
!
! Testing effect of hysteresis 
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
         vpr_diff = e0 - ea(ncell)
!
! Bowen ratio - Andreas_et_al JGR Oceans(2013) Fig 1 
!
         if (T_fit(i) .gt. 0.01) then
           rb = -0.47305*LOG(T_fit(i)) + 1.87629 
         else
           T_rb = T_fit(i) + 40.0
           rb = -5.85894*LOG(T_rb) - 23.3993
         end if
!
! Latent heat of vaporization - joules (Wsec)/kg
         lvp = kcal_Wsec * (597.0 - T_fit(i))
!
! Evaporative heat flux
         QEVAP=wind_fctr*rho*lvp*evap_rate*WIND(ncell)
         if(QEVAP.lt.0.0) QEVAP=0.0
!
! Convective heat flux
         QCONV=rb*QEVAP
         QEVAP=QEVAP*vpr_diff
! 
! Back radiation from the water surface - W/m**2
         QWS=280.23+6.1589*T_fit(i)
!
! Thermal energy budget for i = 1,2
         q_fit(i)=QNS(ncell)+0.97*QNA(ncell)-QWS-QEVAP+QCONV
!if (nd .le.5) write(*,*) 'Heat ',i,nd,ncell,T_fit(i),vpr_diff,lvp &
!                                 ,rb
      end do
!
! Linear relationship for equilibrium temperature and rate constant
!
      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))         &
          /(T_fit(1)-T_fit(2))
!
      qsurf=0.5*(q_fit(1)+q_fit(2))
      RETURN
      END