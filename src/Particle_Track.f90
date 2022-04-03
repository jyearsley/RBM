SUBROUTINE Particle_Track(nr,ns,nx_s)
use Block_Energy
USE Block_Hydro
USE Block_Network
IMPLICIT NONE
integer                     :: nr,ns
integer                     :: jtrp1,jtrp2
integer                     :: ncell, nx_part
integer                     :: nss,nssdmm,nx_s,count_step
logical                     :: DONE_PART
!
real                        :: dt_left,dt_total,dt_dummy,x_dummy
real,dimension(ns_max)      :: dt_sum,xpprt

!
!     Segment is in cell SEGMENT_CELL(NC)
!
!
                    write(*,*) 'Starting here'
                    ncell = segment_cell(nr,ns)
                    if (dt(ncell) .gt. dt_comp) then
                      nstrt_elm(ns) = ns
                      x_part(ns) = x_dist(nr,ns) + u(ncell)*dt_comp
                      nx_s = nx_s + 1
                      jtrp1 = ns
                      jtrp2 = ns-1
                      x1 = x_dist(nr,jtrp1)
                      x2 = x_dist(nr,jtrp2)
                      t1 = temp(nr,jtrp1,n1)
                      t2 = temp(nr,jtrp2,n1)
                       no_dt(ns) = 1
                      dt_part(nx_s) = dt_comp
                      DONE_PART = .TRUE.
!if (nr.eq.1) write(26,*) '0',nr,ns,x_part(ns),x1,x2
                    end if 
!
!                     
 100                continue
       
!     Determine if the total elapsed travel time is equal to the
!     computational interval
!
<<<<<<< HEAD
                    if (.not. DONE_PART) then
                      dt_dummy = 0.0
                      dt_total = 0.0
                      x_dummy = x_dist(nr,ns)
                      do nss = ns,1,-1
                        nssdmm = nss
                        ncell = segment_cell(nr,nss)
                        dt_dummy = dt_dummy + dt(ncell)
                        dt_sum(nss) = dt_dummy
                        dt_total = dt(ncell) + dt_total
                        nx_s = nx_s +  1
                        dt_part(nx_s) = dt(ncell) 
                        nstrt_elm(ns) = nss
                        x_dummy = x_dummy + u(ncell)                 &                      
                                * dt(ncell)
!if (nr.eq.1) write(26,*) '2',ns,nss,dt_part(nx_s) 
                       xpprt(nssdmm) = x_dummy
                        if (dt_total .gt. dt_comp)  then
                          exit
                        end if

                       end do
=======
    do while(dt_total.lt.dt_comp.and.nx_part.gt.0)
       nx_s=nx_s+1
       x_part(nx_s)=x_dist(nr,nx_part-1)
       ncell=segment_cell(nr,nx_part)       
      if (ncell .eq. 0) write(*,*) 'ncell - 0 ',nr,nx_part,ncell
       dt_part(nx_s)=dt(ncell)
       dt_total=dt_total+dt_part(nx_s)
>>>>>>> origin
!
                       if (nstrt_elm(ns) .eq.1 .and.                 &
                          dt_total.lt.dt_comp) then          
                         x_part(ns) = x_dist(nr,0) 
                          DONE_PART = .TRUE.
                       end if                 
!  
                       dt_left = dt_comp - dt_sum(nssdmm+1) 
                       if (.not. DONE_PART .and. dt_left .gt. 0.0) then
                         dt_part(nx_s) = dt_left                 
                         x_part(ns) = xpprt(nssdmm+1)+u(ncell)*dt_left
                       end if
                      if (x_part(ns) .gt. x_bndry) then
                         x_part(ns)=x_dist(nr,0)     
                      end if   
!   
                      jtrp1 = nssdmm
                      jtrp2 = nssdmm-1
                      x1 = x_dist(nr,jtrp1)
                      x2 = x_dist(nr,jtrp2)
                      t1 = temp(nr,jtrp1,n1)
                      t2 = temp(nr,jtrp2,n1)
!
!                      
                    end if
                      DONE_PART = .TRUE.
!if (nr.eq.1) write(26,*) '*****************************'
!
!     For the last segment of particle travel, adjust the particle location
!     such that the total particle travel time is equal to the computational
!     interval.
!   
                    no_dt(ns)=nx_s
!
END SUBROUTINE Particle_Track
!END MODULE P_Track
