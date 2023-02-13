      SUBROUTINE ENERGY_ICE   (q_rslt, nd, ncell, nr, ns)
!
Use Block_Energy
Use Block_Hydro
Use Block_Network
Use Block_Ice_Snow
!
Implicit None
integer                   :: ncell,nd,nr,ns
!
real                      :: cnddmmy,cnd_tmp
real                      :: q_ice,q10
real                      :: q_rslt
real                      :: Ltnt_Heat,Sens_Heat
real                      :: cndctvy,delta_ice,LW_back,LW_in,SW_in
real                      :: dvsr,delta_Temp
real                      :: T_B,T_ice,T_p,T_p_cubed,T_riv,T_srfc
!
real,parameter            :: thick0 = 0.8
!
  T_B = 0.0
  T_srfc = ice_temp(nr,ns,n1)
  T_p = T_srfc + T_Kelvin
  T_p_cubed = T_p*T_p*T_p
!  cndctvy  = ice_cndctvy/(ice_thick(nr,ns,n1)+0.01)
  cndctvy = ice_cndctvy/thick0

  dvsr = 4.0*epsilon*Stfn_Bltz*T_p_cubed + cndctvy 
  LW_back = epsilon*Stfn_Bltz*T_p*T_p_cubed 
  LW_in = epsilon*QNA(ncell)
  SW_in = (1.0-0.4*I0)*(1.0-ice_albedo)*QNS(ncell)

!
  call Vapor(q_ice,q10,T_p,ncell,ns)
  Ltnt_Heat = rho_air*lvs*Ch_Cg*wind(ncell)*(q10-q_ice)
!
  Sens_Heat = rho_air*cp_air*Ch_Cg*wind(ncell)*(dbt(ncell)-T_srfc)
!
  delta_Temp = (Sens_Heat + Ltnt_Heat + LW_in + SW_in - LW_back           &
            +cndctvy*(T_B-T_srfc))/dvsr
  cnd_tmp = cndctvy*(T_B-T_srfc)


  ice_temp(nr,ns,n2) = ice_temp(nr,ns,n1) + delta_Temp
  ice_temp(nr,ns,n2) = AMAX1(-.01,ice_temp(nr,ns,n2))
!
!
! Check here to see if ice is thawing
!
  if (ice_temp(nr,ns,n2) .ge. T_B) then
    
    delta_ice = dt_comp*(Sens_Heat + Ltnt_Heat + LW_in + SW_in - LW_back)/lvf
   ice_thick(nr,ns,n2) = ice_thick(nr,ns,n1) - delta_ice                
    ice_thick(nr,ns,n2) = AMAX1(ice_thick(nr,ns,n2),0.0)
    ice_temp(nr,ns,n2) = AMIN1(ice_temp(nr,ns,n2),T_B)
    ice_temp(nr,1:ns,n2) = T_B
    ICE(ncell) = 102.
!
! Freezing
!
  else
    T_ice = ice_temp(nr,ns,n1)
    T_riv = temp(nr,ns,n1)
    delta_ice = dt_comp*cndctvy*(T_riv-T_ice)/lvf

    ice_thick(nr,ns,n2) = ice_thick(nr,ns,n1) + delta_ice
    ice_thick(nr,ns,n2) = AMIN1(ice_thick(nr,ns,n2),depth(ncell))
    ICE(ncell) = 200.
  end if
!
! If ice thickness is less than the minimum, river is no longer frozen
!
    if (ice_thick(nr,ns,n2) .lt. 0.01) then 
      ICE(ncell) = 101.
      ice_thick(nr,ns,n2) = 0.00
    end if
!
!
  if (ice_thick(nr,ns,n2) .gt. depth(ncell)) ice_thick(nr,ns,n2) = depth(ncell)
!
  q_rslt =  SW_in + LW_in -LW_back - Ltnt_Heat + Sens_Heat

  
      RETURN
      END
