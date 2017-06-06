SUBROUTINE SYSTMM
!
use Block_Energy
use Block_Hydro
use Block_Network
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
real         :: Cl_0,Cl_dist,T_0,T_dist,Cl_00
real(8)      :: time
real         :: x,x_bndry,xd,xdd,xd_year,x_head,xwpd
!
! Indices for lagrangian interpolation
!
integer:: npndx,ntrp

integer, dimension(2):: ndltp=(/-1,-2/)
integer, dimension(2):: nterp=(/2,3/)
!
!
real             :: tntrp
real,dimension(4):: cla,ta,xa
!
real             :: Cl_dist_load,Cl_point_load,Cl_trib_load
real             :: T_dist_load,T_point_load,T_trib_load
real             :: Q_in_mps,Q_out_mps,Q_dist_mps,Q_trib_mps,Q_trib_sum,Q_ratio
real             :: Cl_loading  !Temporary variable for testing
real             :: QQ_in_mps
!
! Allocate chloride
!
allocate (tds(nreach,-2:ns_max,2))
chlr(:,:,:) = 0.0
allocate (tds_trib(nreach))
!
!  Allocate water temperature
!
allocate (temp(nreach,-2:ns_max,2))
temp(:,:,:) = 4.0
allocate (T_head(nreach))
allocate (T_smth(nreach))
allocate (T_trib(nreach))
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
T_head=4.0
!!
!
! Initialize smoothed air temperatures for estimating headwaters temperatures
!
T_smth=4.0
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
        T_smth(nr)=rminsmooth*T_smth(nr)+smooth_param(nr)*dbt(nc_head)
!     
!     Variable Mohseni parameters (UW_JRY_2011/06/16)
! 
        T_head(nr)=mu(nr)+(alphaMu(nr) &
                  /(1.+exp(gmma(nr)*(beta(nr)-T_smth(nr)))))  
!
!  Initialize headwaters temperatures
!
      temp(nr,0,n1)=T_head(nr)
      temp(nr,1,n1)=T_head(nr)
!
      Cl_head(nr) = 0.0
!  Initialize headwaters chlorides
!
      chlr(nr,0,n1) = Cl_head(nr)
      chlr(nr,1,n1) = Cl_head(nr)
!
        do ns=1,no_celm(nr)
! 
          TRIBS_DONE=.FALSE.
          SRCES_DONE=.FALSE.
  
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
          if(nx_head .eq. 0) then
            Cl_0 = Cl_head(nr)
            T_0  = T_head(nr)
          else 
!
!     Interpolation at the upstream or downstream boundary
!
            if(nseg .eq. 1 .or. nseg .eq. no_celm(nr)) npndx=1
!
            do ntrp=nterp(npndx),1,-1
              npart=nseg+ntrp+ndltp(npndx)
              xa(ntrp)=x_dist(nr,npart)
!
!  Interpolate chloride
!            
              cla(ntrp) = chlr(nr,npart,n1)
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
            T_0  = tntrp(xa,ta,x,nterp(npndx))
!            Cl_0 = tntrp(xa,cla,x,nterp(npndx))
            Cl_0 = chlr(nr,nseg,n1)  
            Cl_00 = Cl_0          
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
            z=depth(nncell)
            call energy(T_0,q_surf,nncell)
            q_dot=(q_surf/(z*rfac))
!
!    Add distributed flows
!    
            Cl_dist_load = 0.0
            T_dist_load  = 0.0
            Q_dist_mps = cuft_cum*Q_diff(nncell)
!
            if(Q_diff(nncell).gt.0.001) then
              T_dist  = 10.0
              Cl_dist =  0.0
              T_dist_load  = Q_dist_mps*T_dist
              Cl_dist_load = Q_dist_mps*Cl_dist
            else
              Cl_dist      = Cl_0
              T_dist       = T_0
              Cl_dist_load = Q_dist_mps*Cl_dist
            end if
!
!     Look for a tributary.
!
            ntribs=no_tribs(nncell)
              Q_trib_sum   = 0.0
              Cl_trib_load = 0.0
              T_trib_load = 0.0
!
            if(ntribs.gt.0.and..not.TRIBS_DONE) then
              do ntrb=1,ntribs
                nr_trib=trib(nncell,ntrb)
                if(Q_trib(nr_trib).gt.0.0) then
                  Q_trib_mps = cuft_cum*Q_trib(nr_trib)
                  Q_trib_sum = Q_trib_sum + Q_trib_mps
! 
!  Update chloride with tributary input
!                 
                  Cl_trib_load = (Q_trib_mps*Cl_trib(nr_trib))    &
                               + Cl_trib_load
! 
!  Update water temperature with tributary input
!                 
                  T_trib_load = (Q_trib_mps*T_trib(nr_trib))       &
                              +  T_trib_load
                end if
!
              end do
              TRIBS_DONE = .TRUE.
            end if
!
!  Update inflow and outflow

            Q_out_mps = Q_in_mps + cuft_cum*Q_diff(nncell)
            Q_out_mps = Q_out_mps + Q_trib_sum
            Q_ratio = Q_in_mps/Q_out_mps       
!
            Cl_point_load = 0.0
            T_point_load  = 0.0
!
            if (.not. SRCES_DONE) then
!
!  Add chloride loading
!
              Cl_point_load = chloride(nncell)
!
              if(Cl_0 .lt. 0.0) Cl_0 = 0.0
!
!  Add thermal loading
!
              T_point_load  = thermal(nncell)
!
              if(T_0 .lt. 0.0) T_0=0.0
              SRCES_DONE = .TRUE.
            end if
!
! Do the mass/energy balance
!
              Cl_loading = Cl_0*Q_in_mps 

            Cl_0 = (Cl_0*Q_in_mps                                            &
                 + Cl_point_load + Cl_dist_load + Cl_trib_load)/Q_out_mps   
         write(29,*) nd,nr,nm,nseg,nncell,Cl_0,Cl_loading    &
                               ,Cl_point_load,Q_in_mps,Q_out_mps
            T_0  = T_0*Q_ratio                                               &
                 + (T_point_load + T_dist_load + T_trib_load)/Q_out_mps      &
                 + q_dot*dt_calc              
            if (T_0.lt.0.5) T_0 =0.5
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
            dt_total=dt_total+dt_calc
            QQ_in_mps = Q_in_mps
            Q_in_mps = Q_out_mps
          end do

!
!  Update chloride and water temperture
!        
            chlr(nr,ns,n2) = Cl_0
!
            temp(nr,ns,n2) = T_0
!
!
!  Update tributary chloride
!
            Cl_trib(nr) = Cl_0            
!
!  Update tributary water temperature
!
	    T_trib(nr)  = T_0
!
!   Write all temperature output UW_JRY_11/08/2013
!   The temperature is output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_11/08/2013)
!
            Cl_loading = Cl_0*Q_out_mps 
!            if (mod(ns,2) .eq. 0) then
              call WRITE(time,nd,nr,ncell,ns,T_0,Cl_0,   &
                        T_head(nr),dbt(ncell),CL_loading,      &
                        QQ_in_mps,Q_out_mps)
!                        Q_in(ncell),Q_out(ncell))
!            end if
!
!     End of computational element loop
!
              end do
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
