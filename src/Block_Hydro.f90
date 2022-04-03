!
! Module for hydraulic characteristics and water quality constituents of the basin
!
module Block_Hydro
    integer, dimension(2000):: no_dt,nstrt_elm
!
    real, dimension(:),   allocatable  :: depth
    real, dimension(:),   allocatable  :: x_part
    real, dimension(:),   allocatable  :: width
    real, dimension(:),   allocatable  :: u
    real, dimension(:),   allocatable  :: dt
    real, dimension(:),   allocatable  :: dx
    real, dimension(:),   allocatable  :: Q_in
    real, dimension(:),   allocatable  :: Q_trib
    real, dimension(:),   allocatable  :: Q_out
    real, dimension(:),   allocatable  :: Q_diff
<<<<<<< HEAD
    real, dimension(:),   allocatable  :: dt_part
    real, dimension(:,:), allocatable  :: Q_nps
=======
    real, dimension(:), allocatable  :: Q_nps
    real, dimension(:), allocatable    :: base_flow
    real, dimension(:), allocatable    :: run_off
>>>>>>> origin
    real, dimension(:,:), allocatable  :: temp_trib
    real, dimension(:,:), allocatable  :: temp_nps,thermal
    real, dimension(:,:), allocatable  :: x_dist
    real, dimension(:,:,:), allocatable :: temp

end module Block_Hydro
