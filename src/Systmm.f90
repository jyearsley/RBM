SUBROUTINE SYSTMM
!
use Block_Energy
use Block_Hydro
use Block_Ice_Snow
use Block_Network
!
Implicit None
! 
!
integer          :: ncell, nc_head
integer          :: nd, ndd, nr, ns
integer          :: nobs, nyear, nd_year, ntmp
!
real             :: hpd, xd, xdd, xwpd, xd_year, year
real             :: rminsmooth
real(8)          :: time
!
! Allocate the arrays
!
!allocate (ICE(heat_cells))
!allocate (SNOW(heat_cells))
!allocate (SUB_ZERO(heat_cells))
allocate (T_head(nreach))
allocate (T_smth(nreach))
allocate (T_trib(nreach))
allocate (depth(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_out(heat_cells))
allocate (Q_diff(heat_cells))
allocate (Q_trib(nreach))
allocate (width(heat_cells))
allocate (u(heat_cells))
allocate (dt(2*heat_cells))
allocate (dbt(heat_cells))
allocate (ea(heat_cells))
allocate (Qns(heat_cells))
allocate (Qna(heat_cells))
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
  if (mod(nyear,4).eq.0) nd_year=366
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
      temp(nr,0,n1)=T_head(nr)
      temp(nr,1,n1)=T_head(nr)
!
! Begin cell computational loop
!
        do ns=1,no_celm(nr)  
!
          ncell = segment_cell(nr,ns)
!
! Check to see if air temperature is above freezing (0.0 deg C)
! 
!    
!write(*,*) 'About to check for ice ',nd,nr,ns,ncell,nc_head
          ICE(ncell)      = .FALSE.
          SUB_ZERO(ncell) = .FALSE.
!
!          if (dbt(ncell) .lt. 0.01) SUB_ZERO(ncell) = .TRUE.
        
!write(*,*) 'Checking for ice',dbt(ncell),ICE(ncell),SUB_ZERO(ncell)
!
          if (.not. SUB_ZERO(ncell) .and. .not.ICE(ncell)) then
!
             call Ice_Free (nd,nr,ns,ncell,nc_head)
          else
             call Frozen (nd,nr,ns,ncell,nc_head)
          end if
!
!   Write all temperature output UW_JRY_11/08/2013
!   The temperature is output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_11/08/2013)
!

        call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),Q_inflow,Q_outflow)
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
end SUBROUTINE SYSTMM
