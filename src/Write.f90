SUBROUTINE WRITE(time,nd,nr,ncell,ns,                 &
                 temp_0,tempload_0,temp_point_load,   &
                 tds_0,tdsload_0,tds_point_load,      &
                 Q_in_mps,Q_out_mps,xkm)
!
Implicit NONE
!
integer:: nd,nr,ncell,ns 
real   :: Cl_0,T_0,T_dist
real(8):: time
real   :: temp_0,tempload_0,tds_0,tdsload_0
real   :: temp_point_load,tds_point_load,Q_in_mps,Q_out_mps,xkm
real   :: cl_source,QQ_in,QQ_out
!
WRITE(25,'(f12.4,4i6,3f10.1,f10.0,f10.1,10f15.1)')                &  
           time,nd,nr,ncell,ns,temp_0,tempload_0,temp_point_load, &
           tds_0,tdsload_0,tds_point_load,                        &
           Q_in_mps,Q_out_mps,xkm
return
end SUBROUTINE WRITE
