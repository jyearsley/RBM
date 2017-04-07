SUBROUTINE WRITE(time,nd,nr,ncell,ns,T_0,Cl_0,T_head,dbt,cl_source,QQ_in,QQ_out)
!
Implicit NONE
!
integer:: nd,nr,ncell,ns 
real   :: Cl_0,T_0,T_dist
real(8):: time
real   :: T_head
real   :: dbt
real   :: cl_source,QQ_in,QQ_out
!
write (25,'(f12.4,4i6,3f8.2,f9.0,f8.1,2f10.1)')                &                                              
            time,nd,nr,ncell,ns,T_0,Cl_0,T_head,cl_source      &
           ,QQ_in,QQ_out
end SUBROUTINE WRITE
