SUBROUTINE ICE_FREE(T_0)

use Block_Energy
use Block_Hydro
use Block_Ice_Snow
!
integer          :: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer          :: nc,nd,ndd,nm,nr,ns
integer          :: nr_trib,ntribs
integer          :: nrec_flow,nrec_heat
integer          :: n1,n2,nnd,nobs,nyear,nd_year,ntmp
integer          :: npart,nseg,nx_s,nx_part,nx_head
use Block_Network
!
Implicit None
!
real                              :: T_0,T_dist
real,dimension(4)                 :: ta,xa 
!
!     Establish particle tracks
!
      call Particle_Track(nr,ns,nx_s,nx_head)
!
          ncell=segment_cell(nr,ns)
!
!     Now do the third-order interpolation to
!     establish the starting temperature values
!     for each parcel
!
          nseg=nstrt_elm(ns)
!
!     Perform polynomial interpolation
!
!
!     Interpolation inside the domain
!
          npndx=2
!
!     Interpolation at the upstream boundary if the
!     parcel has reached that boundary
!
          if(nx_head.eq.0) then
            T_0 = T_head(nr)
          else 
!
!
!     Interpolation at the upstream or downstream boundary
!
            if(nseg .eq. 1 .or. nseg .eq. no_celm(nr)) npndx=1
!
            do ntrp=nterp(npndx),1,-1
              npart=nseg+ntrp+ndltp(npndx)
              xa(ntrp)=x_dist(nr,npart)
              ta(ntrp)=temp(nr,npart,n1)
            end do
!
! Start the cell counter for nx_s
!
            x=x_part(nx_s)
!
!     Call the interpolation function
!
            T_0=tntrp(xa,ta,x,nterp(npndx))
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
          ncell0=nncell
          dt_total=0.0
          do nm=no_dt(ns),1,-1
            dt_calc=dt_part(nm)
            z=depth(nncell)
            call energy(T_0,q_surf,nncell)
!
            q_dot=(q_surf/(z*rfac))
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
!            if(ntribs.gt.0.and..not.DONE) then
!
! Uses first segment of the cell to advect tributary thermal energy
!
            min_seg = first_seg(nncell)
            if(ntribs.gt.0.and.nseg.eq.min_seg) then
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
            Q_outflow = Q_inflow + Q_dstrb + Q_trb_sum
            Q_ratio = Q_inflow/Q_outflow       
!
! Do the mass/energy balance
!
            T_0  = T_0*Q_ratio                              &
                 + (T_dstrb_load + T_trb_load)/Q_outflow    &
                 + q_dot*dt_calc              
!
            if (T_0.lt.0.5) T_0 =0.5
            Q_inflow = Q_outflow
!
            nseg=nseg+1
            nncell=segment_cell(nr,nseg)
!
!     Reset tributary flag if this is a new cell
!
            if(ncell0.ne.nncell) then
              ncell0=nncell
              Q_inflow = Q_in(nncell)
               DONE=.FALSE.
            end if
            dt_total=dt_total+dt_calc
          end do
          if (T_0.lt.0.5) T_0=0.5
          temp(nr,ns,n2)=T_0
	      T_trib(nr)=T_0
!
       
END SUBROUTINE ICE_FREE