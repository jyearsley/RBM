SUBROUTINE WRITE(time,nd,nr,ncell,ns,temp_0,tds_0,    &
                        T_head,dbt,tds_loading,       &
                        Q_in_mps,Q_out_mps,xkm)
!
Implicit NONE
!
integer:: nd,nr,ncell,ns 
real   :: Cl_0,T_0,T_dist
real(8):: time
real   :: temp_0,tds_0,T_head,dbt,tds_loading,Q_in_mps,Q_out_mps,xkm
real   :: cl_source,QQ_in,QQ_out
!
WRITE(25,'(f12.4,4i6,3f10.1,f10.0,f10.1,5f15.1)')     & 
           time,nd,nr,ncell,ns,temp_0,tds_0,       &
           T_head,dbt,tds_loading,                 &
           Q_in_mps,Q_out_mps,xkm
return
end SUBROUTINE WRITE
