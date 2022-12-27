SUBROUTINE WRITE(time,nd,nr,ncell,ns,T_0,T_head,dbt,ice_thick,ice_temp,ICC)
!
Implicit NONE
!
integer :: nd,nr,ncell,ns 
logical :: ICC
real    :: Q_inflow,Q_outflow
real    :: ice_thick,ice_temp
real    :: T_0,T_head
real(8) :: time
real    :: dbt
!
write (20,'(f12.4,4i6,3f8.2,f8.4,f8.2,2x,l8)')           &                                              
            time,nd,nr,ncell,ns,T_0,T_head,dbt,ice_thick,ice_temp,ICC
end SUBROUTINE WRITE
