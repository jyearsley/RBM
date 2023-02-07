SUBROUTINE ICE_FREE(nd,nr,ns,ncell,nx_head)

use Block_Energy
use Block_Hydro
use Block_Ice_Snow
use Block_Network
!
Implicit None
!
integer          :: nd,ncell,nncell,ncell0,strt_cell
integer          :: nm,nr,ns
integer          :: nr_trib,ntribs
integer          :: nseg,nx_s,nx_head
!
logical          :: DONE
!
real             :: dt_calc,dt_total,q_dot,q_surf,z
real             :: Q_dstrb,Q_ratio,Q_trb,Q_trb_sum
real             :: T_dstrb,T_dstrb_load,T_trb_load,t_total
real,parameter   :: ft_to_m =1./3.2808
!
DONE = .FALSE.
!
!     Establish particle tracks
!
nx_s = 0
!  
      call Particle_Track(nr,ns,nx_s)
!
          ncell=segment_cell(nr,ns)
          nseg=nstrt_elm(ns)
!
! Check if particle is at x_bndry
!
          if(x_part(ns).gt.x_bndry) then
            T_0 = T_head(nr)
          else 
!
! Linear interpolation
!
            dlta1 = x2 - x_part(ns)
            dlta2 = x_part(ns) - x1
            dlta  = x2 - x1
            T_0 = (dlta1*t1 + dlta2*t2)/dlta
!
          end if
!
!
          nncell=segment_cell(nr,nstrt_elm(ns))
!

!    Initialize inflow
!
          Q_inflow = Q_in(nncell)
          Q_outflow = Q_out(nncell)
!
!    Set NCELL0 for purposes of tributary input
!
          t_total = 0.0
          ncell0=nncell
          dt_total=0.0
          do nm=no_dt(ns),1,-1
            dt_calc=dt_part(nm)
            t_total = t_total + dt_calc
            z=depth(nncell)*ft_to_m
            z = AMAX1(z,2.0)
     if (ncell .eq. 1040 .and. ns .eq. 122) write(86,*)  'Track 2 ',nd,ncell,ncell0 &
                                                         ,nseg,T_0,dt_calc,t_total,z
            call Energy_NO_ICE(T_0,q_surf,nd,nr,nncell,nseg)
!
            q_dot=(q_surf/(z*rho_Cp))
            T_0=T_0+q_dot*dt_calc
            if(T_0.lt.0.0) T_0=0.0
!
!    Add distributed flows
!    
            T_dstrb_load  = 0.0
!
            Q_dstrb = Q_diff(nncell)
!

! Temperature of distributed inflow assumed = 10.0 deg C
!
            if(Q_dstrb.gt.0.001) then
              T_dstrb  = 10.0
            else
              T_dstrb  = 10.0
            end if
              T_dstrb_load  = Q_dstrb*T_dstrb

!
!     Look for a tributary.
!
            ntribs=no_tribs(nncell)
            Q_trb_sum   = 0.0
            T_trb_load  = 0.0
!
! Uses first segment of the cell to advect tributary thermal energy
!
            if(ntribs.gt.0.and..not.DONE) then
!
              do ntrb=1,ntribs
                nr_trib=trib(nncell,ntrb)
                if(Q_trib(nr_trib).gt.0.0) then
                  Q_trb        = Q_trib(nr_trib)
                  Q_trb_sum    = Q_trb_sum + Q_trb
!
!  Update water temperature with tributary input
!
                  T_trb_load   = (Q_trb*T_trib(nr_trib))       &
                               +  T_trb_load
                end if
              end do
!
              DONE=.TRUE.
            end if
!
!  Update inflow and outflow
!
            Q_outflow = Q_inflow + Q_dstrb + Q_trb_sum + 0.0001
            Q_ratio = Q_inflow/Q_outflow       
!
! Do the mass/energy balance
!
            T_0  = T_0*Q_ratio                              &
                 + (T_dstrb_load + T_trb_load)/Q_outflow    &
                 + q_dot*dt_calc                           
!
            Q_inflow = Q_outflow
!
            nseg=nseg+1
            nncell=segment_cell(nr,nseg)
!
!     Reset tributary flag if this is a new cell
!
            if(ncell0.ne.nncell .and. nncell .gt.0) then
              ncell0=nncell
              Q_inflow = Q_in(nncell)
               DONE=.FALSE.
            end if
            dt_total=dt_total+dt_calc
          end do
           if (T_0.lt.0.001) then
            T_0 =0.01
            strt_cell = segment_cell(nr,1) 
            ice_thick(nr,1:ns,n2) = 0.01
            ice_temp(nr,1:ns,n2) = -0.01
          end if
          temp(nr,ns,n2)=T_0
          T_trib(nr)=T_0
!
       
END SUBROUTINE ICE_FREE
