SUBROUTINE WRITE(time,nd,nr,ncell,ns,T_0,T_head,dbt,dpth,ice_thick,ice_temp,ICC)
!
Implicit NONE
!
integer :: nd,nr,ncell,ns 
real    :: ICC
real    :: dpth,Q_inflow,Q_outflow
real    :: ice_thick,ice_temp
real    :: T_0,T_head
real(8) :: time
real    :: dbt
!
write (20,'(f12.4,4i6,4f8.2,f8.1,f8.2,2x,f8.1)')           &                                              
            time,nd,nr,ncell,ns,T_0,T_head,dbt,dpth,ice_thick,ice_temp,ICC
end SUBROUTINE WRITE
