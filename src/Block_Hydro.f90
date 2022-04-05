!
! Module for hydraulic characteristics and water quality constituents of the basin
!
module Block_Hydro
!
    integer, dimension(2000):: no_dt,nstrt_elm
!
!
    real                    :: Q_inflow, Q_outflow
    real, dimension(2000)   :: dt_part,x_part
!
    real, dimension(:),   allocatable  :: depth
    real, dimension(:),   allocatable  :: width
    real, dimension(:),   allocatable  :: u
    real, dimension(:),   allocatable  :: dt
    real, dimension(:),   allocatable  :: dx
    real, dimension(:),   allocatable  :: Q_in
    real, dimension(:),   allocatable  :: Q_trib
    real, dimension(:),   allocatable  :: Q_out
    real, dimension(:),   allocatable  :: base_flow
    real, dimension(:),   allocatable  :: run_off
    real, dimension(:),   allocatable  :: Q_diff

    real, dimension(:,:), allocatable  :: x_dist


end module Block_Hydro
