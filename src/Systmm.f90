SUBROUTINE SYSTMM(temp_file,param_file)
!
use Block_Energy
use Block_Hydro
use Block_Network
!
Implicit None
! 
!
character (len=200):: temp_file
character (len=200):: param_file
! 
integer          :: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer          :: nc,nd,ndd,nm,nr,ns
integer          :: nr_trib,ntribs
integer          :: nrec_flow,nrec_heat
<<<<<<< HEAD
integer          :: nnd,nobs,nyear,nd_year,ntmp
integer          :: npart,nseg,nx_s,nx_part
=======
integer          :: n1,n2,nnd,nobs,nyear,nd_year,ntmp,nrr_tmp
integer          :: npart,nseg,nx_s,nx_part,nx_head
>>>>>>> origin
!
! Indices for lagrangian interpolation
!
integer              :: njb,npndx,ntrp
integer, dimension(2):: ndltp=(/-1,-2/)
integer, dimension(2):: nterp=(/2,3/)

!
real             :: dt_calc,dt_total,hpd,q_dot,q_surf,z
real             :: Q_dstrb,Q_inflow,Q_outflow,Q_ratio,Q_trb,Q_trb_sum

real             :: T_dstrb,T_dstrb_load,T_trb_load
real             :: rminsmooth
real             :: Qload1,T_0,T_dist
real(8)          :: time
real             :: x,xd,xdd,xd_year,xwpd,year
real             :: tntrp
real             :: dt_ttotal
real,dimension(4):: ta,xa
!
real,dimension(:),allocatable     :: T_head,T_smth,T_trib

logical:: DONE
!
logical Leap_Year
!
! Allocate the arrays
!
allocate (dt_part(ns_max))
allocate (x_part(ns_max))
allocate (T_head(nreach))
allocate (T_smth(nreach))
allocate (T_trib(nreach))
allocate (depth(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_out(heat_cells))
allocate (Q_diff(heat_cells))
allocate (base_flow(heat_cells))
allocate (run_off(heat_cells))
allocate (Q_trib(nreach))
allocate (width(heat_cells))
allocate (u(heat_cells))
allocate (dt(2*heat_cells))
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
!
!
! Initialize smoothed air temperatures for estimating headwaters temperatures
!
T_smth=4.0

!
!     open the output file
!

open(20,file=TRIM(temp_file),status='unknown')
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
!  if (mod(nyear,4).eq.0) nd_year=366
  if (Leap_Year(nyear)) nd_year = 366
!
!     Day loop starts
!
  DO nd=1,nd_year
    year=nyear
    xd=nd
    xd_year=nd_year
!     Start the numbers of days-to-date counter
!
    ndays=ndays+1
!
!    Daily period loop starts
!
      DO ndd=1,nwpd
      xdd = ndd
<<<<<<< HEAD
      time=year+(xd+(xdd-0.5)*hpd)/xd_year 
=======
      time=year+(xd+(xdd-1.5)*hpd)/xd_year 
>>>>>>> origin

!
! Read the hydrologic and meteorologic forcings
!
        call READ_FORCING
        nrr_tmp = 0
        do nr = 1,nreach
            do nc=1,no_cells(nr)
            nrr_tmp = nrr_tmp + 1
            end do
        end do

!
!     Begin reach computations
!
!
!     Begin cycling through the reaches
!
      do nr=1,nreach
!
      nrr_tmp = nr
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
      temp(nr,0,n1)=T_head(nr)
      temp(nr,1,n1)=T_head(nr)
      x_bndry=x_dist(nr,0) - 1.0      
!
      DONE=.FALSE.
!
      DONE=.FALSE.
!
! Begin cell computational loop
!
        do ns=1,no_celm(nr)
! 
! Testing new code 8/8/2016
!
!     Establish particle tracks
!
      nx_s = 0

!
<<<<<<< HEAD
      call Particle_Track(nr,ns,nx_s)
      ncell=segment_cell(nr,ns)  
=======
          ncell=segment_cell(nr,ns)
!
!     Now do the third-order interpolation to
!     establish the starting temperature values
!     for each parcel
>>>>>>> origin
!
! Start the cell counter for nseg
!          
          nseg=nstrt_elm(ns)
!
! Check if particle is at x_bndry
!
          if(x_part(ns).gt.x_bndry) then
            T_0 = T_head(nr)
!write(26,*) 'Headwaters',nr,nseg,T_0
          else 
!
! Do the interpolation
!
<<<<<<< HEAD
            dlta1 = x2 - x_part(ns)
            dlta2 = x_part(ns) - x1
            dlta = x2 - x1           
            x=x_part(ns)
            T_0 = (dlta1*t1 + dlta2*t2)/dlta
!if (nr.eq.1) write(26 ,*) 'Interp ',T_0,t1,t2

=======
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
>>>>>>> origin
          end if
!
!
          nncell=segment_cell(nr,nstrt_elm(ns))
!
!    Set NCELL0
!
          ncell0 = nncell
!
!    Initialize inflow
!
          Q_inflow = Q_in(nncell)
          Q_outflow = Q_out(nncell)
<<<<<<< HEAD
!
!    Set NCELL0 for purposes of tributary input
!
=======
!          Q_outflow = Q_in(nncell) +Q_diff(nncell)
          
>>>>>>> origin
          dt_total=0.0
          do nm=no_dt(ns),1,-1
            dt_calc=dt_part(nm)
            z=depth(nncell)
<<<<<<< HEAD
!if (nr.eq.1) write(26,*) 'Heat',nr,ns,nm,nncell,T_0,z
            call energy(T_0,q_surf,nncell,nr)
!
            q_dot=(q_surf/(z*rfac))
!
! The following update for T_0 is redundant per RJN 7/26/3017
=======
!
            call energy(T_0,q_surf,nncell,nrr_tmp)
!
            q_dot=(q_surf/(z*rfac))
!
! The following update for T_0 is redundant per RJN 7/26/2017
>>>>>>> origin
! and has been commented out for the time being - JRY
!
!            T_0=T_0+q_dot*dt_calc
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
            Q_outflow = Q_inflow + Q_dstrb + Q_trb_sum
            Q_ratio = Q_inflow/(Q_outflow + 0.1)          ! Temporary fix JRY 06/24/2021      
!
! Do the mass/energy balance
!
            T_0  = T_0*Q_ratio                              &
                 + (T_dstrb_load + T_trb_load)/(Q_outflow + 0.1)    &
                 + q_dot*dt_calc  
                 
   Qload1 = q_dot*dt_calc

!
            if (T_0.lt.0.5) T_0 =0.5
            Q_inflow = Q_outflow
!
<<<<<<< HEAD
            ncell0=nncell
=======
!            ncell0=nncell
>>>>>>> origin
            nseg=nseg+1
            nncell=segment_cell(nr,nseg)
!
!     Reset tributary flag if this is a new cell
!
            if(ncell0.ne.nncell) then
              ncell0=max(nncell,ncell0)
              Q_inflow = Q_in(ncell0)
               DONE=.FALSE.
            end if
            dt_total=dt_total+dt_calc
          end do
          if (T_0.lt.0.5) T_0=0.5
            temp(nr,ns,n2)=T_0
	    T_trib(nr)=T_0
!
!   Write all temperature output UW_JRY_11/08/2013
!   The temperature is output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_11/08/2013)
!
<<<<<<< HEAD
            call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),Q_inflow,Q_outflow)
=======
  call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell) &
                      ,Q_inflow,Q_outflow)
>>>>>>> origin
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
