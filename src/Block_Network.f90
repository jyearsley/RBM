Module Block_Network
!
! Module with stream topology variables
!
    integer, dimension(:), allocatable  :: no_celm,no_cells,no_tribs
    integer, dimension(:), allocatable  :: head_cell
!
    integer, dimension(:,:), allocatable::conflnce,reach_cell,segment_cell,trib
!
!
! Integer variables 
!
    integer           :: flow_cells,heat_cells
    integer           :: ndays,nreach,ntrb,nwpd
    integer,parameter :: ns_max=500
    integer           :: start_year,start_month,start_day
    integer           :: end_year,end_month,end_day
!
! Real variables
!
    real                           :: delta_n,dt_comp
    real                           :: year
    real                           :: n_default = 2.0
    real, dimension(:),allocatable :: ndelta
end module Block_Network