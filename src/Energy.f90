SUBROUTINE Energy(T_surf,q_surf,ncell,nr)
   use Block_Energy
   use Block_Network
   implicit none
   integer           :: i,ncell,nd,nr
   real              :: A,B,e0,q_surf,vpr_diff
   real              :: HV_in,LW_back,LV_in,LW_in,SW_in
   real              :: td,T_rb,T_surf,T_Kelvin,T_tetens
   real, dimension(2):: q_fit, T_fit
!
   td=nd
   T_fit(1)=T_surf-1.0
   T_fit(2)=T_surf+1.0
   do i=1,2
      T_Kelvin = T_fit(i) + 273.0
      T_tetens = T_fit(i) + 237.3
!
! Vapor pressure at water surface -kPa (Magnus-Tetens)
!
      e0=0.61078*exp((17.27*T_fit(i))/T_tetens)
!
!Bowen ratio - Andreas_et_al JGR Oceans (2013) Figure 1
!
      T_rb = DBT(ncell)
      if (T_rb .lt. 0.05) T_rb = 0.05
!
      lvp = kcal_Wsec*(597.0-T_fit(i))
      vpr_diff = e0 - ea(ncell)
!
! Evaporative head flux - uses only the Lake Hefner coefficient
      LV_in = wind_fctr*rho_H2O*lvp*evrate(1)*wind(ncell)
      LV_in = LV_in*vpr_diff
!      if(LV_in .lt. 0.0) LV_in=0.0
!
! Convective heat flux
      HV_in=rb*(a_conv*wind(ncell)+b_conv)*(dbt(ncell)-T_fit(i))
!
! Shortwave radiation
      SW_in = Q_NS(ncell)
!
! Longwave radiation
      LW_in = Q_NA(ncell)
!
! Back radiation
      LW_back = 280.23 + 6.1589*T_fit(i)
!
      q_fit(i) = SW_in + LW_in - LW_back - LV_in + HV_in    
!
   end do
!
!     q=AT+B
!
!     Linear fit over the range of 2.0 deg C.
!     These results can be used to estimate the "equilibrium" 
!     temperature and linear rate constant.
!
   A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
   q_surf=0.5*(q_fit(1)+q_fit(2))
   B=(q_surf/A)-(T_fit(1)+T_fit(2))/2.
!
!     ******************************************************
!               Return to Subroutine RIVMOD
!     ******************************************************
!
END Subroutine Energy
