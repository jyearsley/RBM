SUBROUTINE Particle_Track(nr,ns,nx_s,nx_head,xprt)
USE Block_Hydro
USE Block_Network
IMPLICIT NONE
integer, intent(IN) :: nr,ns
integer             :: ncell, nx_head, nx_part, nx_s
real                :: dtt, dttl, dt_total, dt_xcess
real                :: xprt
!
!     First do the reverse particle tracking
!
!     Segment is in cell SEGMENT_CELL(NC)
   nx_s=0
   nx_part=ns
   dt_xcess = 0.0
   dt_total=0.0
   dtt = 0.0

!
!     Determine if the total elapsed travel time is equal to the
!     computational interval
!
    do while(dt_total.lt.dt_comp.and.nx_part.gt.0)
       nx_s=nx_s+1
!
       ncell = segment_cell(nr,nx_part)       
       xprt = x_dist(nr,nx_part)

!
       if (xprt .gt. x_dist(nr,0)) xprt = x_dist(nr,0)
!
       dt_part(nx_s) = dt(ncell)
!
!
       dttl = dt_total
!       
       dt_total = dt_total+dt_part(nx_s)
       if (dt_total .gt. dt_comp) then
         dt_xcess = dt_total - dt_comp
         dt_part(nx_s) = dt_total - dt_xcess - dttl
         dttl = dttl + dt_part(nx_s)
         xprt = xprt + u(ncell)*dt_part(nx_s)
       end if
!
!     Increment the segment counter if the total time is less than the
!     computational interval
!
        nx_part = nx_part-1
    end do
!
    nx_head=nx_part
    nx_part=max(1,nx_part)
    nstrt_elm(ns)=nx_part
    no_dt(ns)=nx_s 
!
END SUBROUTINE Particle_Track
