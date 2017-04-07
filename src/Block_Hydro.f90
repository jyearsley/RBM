!
! Module for hydraulic characteristics and water quality constituents of the basin
!
module Block_Hydro
    integer, dimension(200):: no_dt,nstrt_elm
    real, dimension(200)   :: dt_part,x_part
!
    real, dimension(:),     allocatable  :: depth
    real, dimension(:),     allocatable  :: width
    real, dimension(:),     allocatable  :: u
    real, dimension(:),     allocatable  :: dt
    real, dimension(:),     allocatable  :: dx
    real, dimension(:),     allocatable  :: Q_in
    real, dimension(:,:),   allocatable  :: Q_in_seg
    real, dimension(:),     allocatable  :: Q_trib
    real, dimension(:),     allocatable  :: Q_out
    real, dimension(:,:),   allocatable  :: Q_out_seg
    real, dimension(:),     allocatable  :: Q_diff
    real, dimension(:,:),   allocatable  :: Q_nps
    real, dimension(:),     allocatable  :: Cl_trib
    real, dimension(:),     allocatable  :: T_smth,T_trib
    real, dimension(:),     allocatable  :: chloride,thermal
    real, dimension(:,:),   allocatable  :: x_dist
    real, dimension(:,:,:), allocatable :: chlr,temp
!
!  Added headwaters variable for chloride here as temporary fix to water quality model
!   
    real, dimension(:),   allocatable  :: T_head,Cl_head

!  Parameters
!
    real,parameter          :: cuft_cum = 0.028318 ! Cubic feet to cubic meters
!
end module Block_Hydro
