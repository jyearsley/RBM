SUBROUTINE WRITE(time,nd,nr,ncell,ns,T_0,T_head,dbt,Q_inflow,Q_outflow,Q_mps)
!
Implicit NONE
!
integer :: nd,nr,ncell,ns 
real    :: Q_inflow,Q_outflow,Q_mps
real    :: T_0,T_dist
real(8) :: time
real    :: T_head
real    :: dbt
!
write (20,'(f12.4,4i6,3f8.2,10f15.1)')           &                                              
            time,nd,nr,ncell,ns,T_0,T_head,dbt,Q_inflow,Q_outflow,Q_mps
end SUBROUTINE WRITE
