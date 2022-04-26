      SUBROUTINE ENERGY   (Tsurf, q_rslt, nd, nr, ncell, ns)
!
   use Block_Energy
   use Block_Ice_Snow
   implicit none
   integer             ::i, ncell, nd, nr, ns
   real                :: A, B, e0, evap_rate,rb,vpr_diff
   real                :: Denom,LV_in,H_in,SW_in,LW_in,LW_back
   real                :: q10,q_ice,q_rslt
   real                :: T_Kelvin,T_p,T_p_cubed,T_rb, Tsurf, T_tetens
   real, dimension(2)  :: q_fit, T_fit
!
! Testing effect of hysteresis 
!
      evap_rate=evrate(1)
      if (nd > 180) evap_rate=evrate(2)
      T_fit(1) = Tsurf-1.0
      T_fit(2) = Tsurf+1.0
!
      do i=1,2
         T_kelvin = T_fit(i) + 273.0
         T_tetens = T_fit(i) + 237.3
!
! Vapor pressure at water surface - kPa (Tetens)
         e0=0.61078*exp((17.27*T_fit(i))/T_tetens)
!
! Bowen ratio - Andreas_et_al JGR Oceans(2013) Figures 1 & 6
!
         if (T_fit(i) .ge. 0.01) then
           T_rb = T_fit(i)
           if (T_rb .lt. 0.01) T_rb = 0.01
          rb = 0.40*(-0.196231*LOG(T_rb)) + 1.411189
!
         else
           T_rb = T_fit(i) + 40.0
           rb = 0.4*(-8.688289*LOG(T_rb)) + 32.232621
         end if
!
! Latent heat of vaporization/sublimation - joules (Wsec)/kg
!         
         if (ICE(ncell)) then
           T_p = ice_temp(nr,ns,1) + 273.0 ! Notation after Parkinson-Washington
           T_p_cubed = T_p*T_p*T_p
! Sensible heat with ice cover
           H_in  = rho_air*cp_air*Ch_Cg*WIND(ncell)*(DBT(ncell)-T_p)           
! Latent heat heat
           q10 = eps_air_h2o*ea(ncell)/(press(ncell)-comp_eps*ea(ncell))
           q_ice = eps_air_h2o*e0/(press(ncell)-comp_eps*e0) 
           LV_in = rho_air*lvs*Ch_Cg*WIND(ncell)*(q10 - q_ice)
! Shortwave radiation
           SW_in = (1.-alpha_ice)*QNS(ncell)
! Longwave radiation 
           LW_in = epsilon*QNA(ncell) 
! Back radiation
           LW_back = epsilon*Stfn_Bltz*T_p*T_p_cubed
! Heat transfer through ice                   
           Ice_Trnsfr = (kappa_ice/ice_thick(nr,ns,1))*(273.0-T_p)
! Denominator
           Denom = 4.0*epsilon*Stfn_Bltz*T_p_cubed + kappa_ice/ice_thick(nr,ns,1)
!
           q_fit(i) = (H_in + LV_in + LW_in + SW_in - LW_back)/Denom
         else 
           lvp = kcal_Wsec * (597.0 - T_fit(i))
           vpr_diff = e0 - ea(ncell)
!
! Evaporative heat flux
         LV_in=wind_fctr*rho*lvp*evap_rate*WIND(ncell)
         if(LV_in.lt.0.0) LV_in=0.0
!
! Convective heat flux
         H_in=rb*LV_in
! Evaporative heat flux
         LV_in=LV_in*vpr_diff
! Shortwave radiation
         SW_in = QNS(ncell)
!Longwave radiaton
         LW_in = QNA(ncell) 
! Back radiation from the water surface - W/m**2
         LW_back=280.23+6.1589*T_fit(i)
!
        end if
!
! Thermal energy budget for i = 1,2
!         q_fit(i)=QNS(ncell)+0.97*QNA(ncell)-QWS-QEVAP+QCONV
         q_fit(i)=SW_in + LW_in - LW_back - LV_in + H_in
         
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
      q_rslt=0.5*(q_fit(1)+q_fit(2))
      RETURN
      END
