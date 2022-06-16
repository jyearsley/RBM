Program Create_Ice
!
Use Block_Energy
Use Block_Hydro
Use Block_Network
Use Block_Ice_Snow
!
Implicit None
!
real                      :: ltnt_heat
real                      :: H,LE
!
ltnt_heat = lvp
if (ICE(ncell)) then
  if (h_rate(nr,ncell)) then
    ltnt_heat = lvf
  else
    ltnt_heat = lvs
  end if
end if
!
H = rho_air*cp_air*Ch_Cg*wind(ncell)*(dbt(ncell)-ice_temp(nr,ncell))
End Program Create_Ice
