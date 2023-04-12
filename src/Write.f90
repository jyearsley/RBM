SUBROUTINE WRITE(time,nd,nr,ncell,ns,T_0,T_head,dbt,dpth,Q_inflow,ice_thick,ICC)
!
Implicit NONE
!
integer        :: nd,nr,ncell,ns 
real           :: ICC
real           :: dpth,Q_inflow,Q_outflow
real           :: ice_thick,ice_temp
real           :: T_0,T_head
real(8)        :: time
real           :: dbt
real,parameter :: cuft_cum = 1.0/(3.32808*3.2808*3.2808)
!
Q_inflow = cuft_cum*Q_inflow
!
write (20,'(f12.4,4i6,4f8.2,f8.1,f8.2,2x,f8.1)')           &                                              
            time,nd,nr,ncell,ns,T_0,T_head,dbt,dpth,Q_inflow,ice_thick,ICC
!
end SUBROUTINE WRITE
