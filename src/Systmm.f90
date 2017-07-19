SUBROUTINE SYSTMM
!
use Block_Energy
use Block_Hydro
use Block_Network
use Block_WQ
!
Implicit None
!
!
logical      :: SRCES_DONE,TRIBS_DONE,LEAP_YEAR
integer      :: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer      :: nc,nd,ndd,nm,nr,ns,ndmmy
integer      :: nr_trib,ntrb,ntribs
integer      :: nrec_flow,nrec_heat
integer      :: n1,n2,nnd,nobs,nyear,nd_year,ntmp
integer      :: npart,nseg,nx_s,nx_head 
real         :: dt_calc,dt_total,hpd,Q1,Q2
real         :: q_dot,q_surf,z
real         :: rminsmooth
real         :: tds_0,tds_dist,temp_0,temp_dist,tds_00,tds_loading
real(8)      :: time
real         :: x,x_bndry,xd,xdd,xd_year,x_head,xwpd
!
! Indices for lagrangian interpolation
!
integer:: npndx,ntrp

integer, dimension(2):: ndltp=(/-1,-2/)
integer, dimension(2):: nterp=(/2,3/)
!
real, dimension(:), allocatable  :: seg_load
!
real             :: tntrp
real, dimension(4):: tdsa,ta,xa
real, dimension(:), allocatable :: temp_smth
!
real             :: tds_dist_load,tds_point_load,tds_trib_load
real             :: temp_dist_load,temp_point_load,temp_trib_load
real             :: Q_in_mps,Q_out_mps,Q_dist_mps,Q_trib_mps,Q_trib_sum,Q_ratio
real             :: Cl_loading  !Temporary variable for testing
real             :: QQ_in_mps
!
! Allocate chloride
!
allocate (tds(nreach,-2:ns_max,2))
tds(:,:,:) = 0.0
allocate (tds_trib(nreach))
allocate (tds_head(nreach))
!
!  Allocate water temperature
!
allocate (temp(nreach,-2:ns_max,2))
temp(:,:,:) = 4.0
allocate (temp_head(nreach))
allocate (temp_smth(nreach))
allocate (temp_trib(nreach))
!
!  Allocate hydrologic input
!
allocate (depth(heat_cells))
allocate (width(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_in_seg(nreach,ns_max))
allocate (Q_out(heat_cells))
allocate (Q_out_seg(nreach,ns_max))
allocate (Q_diff(heat_cells))
allocate (Q_trib(nreach))
allocate (u(heat_cells))
allocate (dt(10*heat_cells))
!
!  Allocate thermal energy budget input
!
allocate (dbt(heat_cells))
allocate (ea(heat_cells))
allocate (Q_ns(heat_cells))
allocate (Q_na(heat_cells))
allocate (press(heat_cells))
allocate (wind(heat_cells))
!
! Initialize some arrays
!
dt_part=0.
x_part=0.
no_dt=0
nstrt_elm=0
temp=0.5
! Initialize headwaters temperatures
!
temp_head=4.0
!!
!
! Initialize smoothed air temperatures for estimating headwaters temperatures
!
temp_smth=4.0
!
!
! Initialize dummy counters that facilitate updating simulated values
!
n1=1
n2=2
nobs=0
ndays=0
xwpd=nwpd
hpd=1./xwpd
!
!     Year loop starts
!
do nyear=start_year,end_year
  write(*,*) ' Simulation Year - ',nyear,start_year,end_year
  nd_year=365
!
! Check to see if it is a leap year
!
  if (LEAP_YEAR(nyear)) nd_year=366
  write(*,*) 'Days in year',nd_year
!
!     Day loop starts
!
  DO nd=1,nd_year
!
    year=nyear
    xd=nd
    xd_year=nd_year
!     
!     Start the numbers of days-to-date counter
!
    ndays=ndays+1
!
!    Daily period loop starts
!
      DO ndd=1,nwpd
      xdd = ndd
      time=year+(xd+(xdd-0.5)*hpd)/xd_year 

!
! Read the hydrologic and meteorologic forcings
!
        call READ_FORCING
!
!     Begin reach computations
!
!
!     Begin cycling through the reaches
!
      do nr=1,nreach

!
        nc_head=segment_cell(nr,1)
!
!     Determine smoothing parameters (UW_JRY_2011/06/21)
!
        rminsmooth=1.0-smooth_param(nr)
        temp_smth(nr)=rminsmooth*temp_smth(nr)         &
                     +smooth_param(nr)*dbt(nc_head)
!     
!     Variable Mohseni parameters (UW_JRY_2011/06/16)
! 
        temp_head(nr)=mu(nr)+(alphaMu(nr) &
                     /(1.+exp(gmma(nr)*(beta(nr)-temp_smth(nr)))))  
!
!  Initialize headwaters temperatures
!
      temp(nr,0,n1)=temp_head(nr)
      temp(nr,1,n1)=temp_head(nr)
      temp_head(nr) = 0.0
!
!  Initialize headwaters chlorides
   
!
      tds_head(nr) = 0.0
      tds(nr,0,n1) = tds_head(nr)
      tds(nr,1,n1) = tds_head(nr)
!
!
        SRCES_DONE = .FALSE.
allocate (seg_load(no_celm(nr)))
seg_load = tds_source(nr,:)
        do ns=1,no_celm(nr)
        write(27,*)
! 
          TRIBS_DONE=.FALSE.

  
! Testing new code 8/8/2016
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
!     Perform polynomial interpolationnr_trib
!
!
!     Interpolation inside the domain
!
          npndx=2
!
!     Use the headwaters value if particle reaches the upstream boundary
!
!          if(nx_head .eq. 0) then
           if (segment_cell(nr,nseg) .eq. head_cell(nr)) then
            tds_0   = tds_head(nr)
            temp_0  = temp_head(nr)
          else 
!
!     Interpolation at the upstream or downstream boundary
!
!            if(nseg .gt. 1 .or. nseg .eq. no_celm(nr)) then
!
            do ntrp=nterp(npndx),1,-1
              npart=nseg+ntrp+ndltp(npndx)
              xa(ntrp)=x_dist(nr,npart)
!
!  Interpolate TDS
!            
!              tdsa(ntrp) = tds(nr,npart,n1)
!
!  Interpolate temperature
!            
              ta(ntrp)=temp(nr,npart,n1)
            end do
!
            x=x_part(nx_s)
!
!     Call the interpolation function
!
!            temp_0  = tntrp(xa,ta,x,nterp(npndx))
!            tds_0 = tntrp(xa,tdsa,x,nterp(npndx))
            temp_0 = temp(nr,nseg,n1)
            tds_0  = tds(nr,ncell-1,n1)
            tds_00 = tds_0          
            end if
300 continue
350 continue
!
          nncell=segment_cell(nr,nstrt_elm(ns))
!
!    Initialize inflow
!
!          Q_in_mps = cuft_cum*Q_in(nncell)
!
! Q_in is the segment flow where the parcel begins
!
          Q_in_mps = cuft_cum*Q_in_seg(nr,nstrt_elm(ns))

!
!    Set NCELL0 for purposes of tributary input
!
          ncell0=nncell
          dt_total=0.0
          do nm=no_dt(ns),1,-1
            dt_calc=dt_part(nm)
            dt_total=dt_total+dt_calc
            z=depth(nncell)
            call energy(temp_0,q_surf,nncell)
            q_dot=(q_surf/(z*rfac))
!
!    Add distributed flows
!    
            tds_dist_load = 0.0
            temp_dist_load  = 0.0
            Q_dist_mps = cuft_cum*Q_diff(nncell)
!
            if(Q_diff(nncell).gt.0.001) then
              temp_dist  = 10.0
              tds_dist =  0.0
              temp_dist_load  = Q_dist_mps*temp_dist
              tds_dist_load = Q_dist_mps*tds_dist
            else
              tds_dist       = tds_0
              temp_dist      = temp_0
              temp_dist_load = Q_dist_mps*temp_dist
              tds_dist_load  = Q_dist_mps*tds_dist
            end if
!
!     Look for a tributary.
!
            ntribs=no_tribs(nncell)
              Q_trib_sum     = 0.0
              temp_trib_load = 0.0
              tds_trib_load  = 0.0
!
!            if(ntribs.gt.0.and..not.TRIBS_DONE) then
            if(ntribs.gt.0.and.nseg.eq.first_seg(nncell)) then
              do ntrb=1,ntribs
                nr_trib=trib(nncell,ntrb)
                if(Q_trib(nr_trib).gt.0.0) then
                  Q_trib_mps = cuft_cum*Q_trib(nr_trib)
                  Q_trib_sum = Q_trib_sum + Q_trib_mps
! 
!  Update chloride with tributary input
!                 
                  tds_trib_load = (Q_trib_mps*tds_trib(nr_trib))    &
                               + tds_trib_load
! 
!  Update water temperature with tributary input
!                 
                  temp_trib_load = (Q_trib_mps*temp_trib(nr_trib))  &
                                 + temp_trib_load
                end if
!
              end do
              TRIBS_DONE = .TRUE.
            end if
!
!  Update inflow and outflow

            Q_out_mps = Q_in_mps + cuft_cum*Q_diff(nncell)
            Q_out_mps = Q_out_mps + Q_trib_sum
            Q_ratio   = Q_in_mps/Q_out_mps       
!
            tds_point_load   = 0.0
            temp_point_load  = 0.0
!
!            if (.not. SRCES_DONE) then
!
!  Add chloride loading
!
!              tds_point_load = tds_source(nr,nseg)
              tds_point_load = seg_load(nseg)
!
              if(tds_0 .lt. 0.0) tds_0 = 0.0
!
!  Add thermal loading
!
              temp_point_load  = temp_source(nr,nseg)
!
              if(temp_0 .lt. 0.0) temp_0=0.0
!              SRCES_DONE = .TRUE.
!            end if
!
! Do the mass/energy balance
!
            tds_00 = tds_0
            tds_0 =  tds_0*Q_ratio                &
                  + (tds_point_load                         &
                  +  tds_dist_load                          &
                  +  tds_trib_load)/Q_out_mps  
                              tds_loading = tds_0*Q_out_mps  
 if(nr.eq.2)          &
 write(27,*) 'Loadings ',ns,nncell,nseg,tds_00,tds_0,tds(nr,nseg,n1) &
                  ,tds_loading,tds_point_load,tds_dist_load,tds_trib_load
!
            temp_0  =  temp_0*Q_ratio                 &
                    + (temp_point_load                &
                    +  temp_dist_load                 &
                    +  temp_trib_load)/Q_out_mps      &
                    +  q_dot*dt_calc              
            if (temp_0.lt.0.5) temp_0 =0.5
500 continue
            nseg=nseg+1
            nncell=segment_cell(nr,nseg)
!
!     Reset tributary and source flags if this is a new cell
!
            if(ncell0.ne.nncell) then
              ncell0=nncell
              TRIBS_DONE = .FALSE.
              SRCES_DONE = .FALSE.
            end if
!
!  Update dt_calc and Q_in_mps
!

            QQ_in_mps = Q_in_mps
            Q_in_mps = Q_out_mps
          end do

!
!  Update chloride and water temperture
!        
            tds(nr,ns,n2)  = tds_0
!
            temp(nr,ns,n2) = temp_0
!
!
!  Update tributary chloride
!
            tds_trib(nr)   = tds_0            
!
!  Update tributary water temperature
!
	    temp_trib(nr)  = temp_0
!
!   Write all temperature output UW_JRY_11/08/2013
!   The temperature is output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_11/08/2013)
!
            tds_loading = tds_0*Q_out_mps 
!            if (mod(ns,2) .eq. 0) then
              call WRITE(time,nd,nr,ncell,ns,temp_0,tds_0,    &
                        temp_head(nr),dbt(ncell),tds_loading, &
                        QQ_in_mps,Q_out_mps)
!                        Q_in(ncell),Q_out(ncell))
!            end if
!
!     End of computational element loop
!
              end do
deallocate (seg_load)
!     End of reach loop
!
            end do
            ntmp=n1
            n1=n2
            n2=ntmp
!
!     End of weather period loop (NDD=1,NWPD)
!
4650 format(16x,12(6x,f6.0,6x))
4700 format(f10.4,f6.0,15(f6.1,f8.3))
4750 format(f10.4,10(i4,f8.0))
          end do
!
!     End of main loop (ND=1,365/366)
!
        end do
!
!     End of year loop
!
      end do
!
!
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
950 return
end SUBROUTINE SYSTMM
