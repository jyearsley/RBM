module block_wq
!
! Dimensioned and allocated water quality variables 
!
      real, dimension(:), allocatable:: DO_head
      real, dimension(:), allocatable:: BOD_head
      real, dimension(:), allocatable:: TEMP_head
      real, dimension(:), allocatable:: PO4_head
      real, dimension(:), allocatable:: P_Org_head
      real, dimension(:), allocatable:: NO2_head
      real, dimension(:), allocatable:: NO3_head
      real, dimension(:), allocatable:: NH4_head
      real, dimension(:), allocatable:: pH_head
      real, dimension(:), allocatable:: H2CO3_head
      real, dimension(:), allocatable:: HCO3_head
      real, dimension(:), allocatable:: CO3_head
      real, dimension(:), allocatable:: ALK_head
      real, dimension(:), allocatable:: CHLR_head
      real, dimension(:), allocatable:: ALGAE_1_head
      real, dimension(:), allocatable:: ALGAE_2_head
      real, dimension(:), allocatable:: ZOO_1_head
      real, dimension(:), allocatable:: ZOO_2_head
      real, dimension(:), allocatable:: TDS_head

!
! Dimensioned and allocated water quality variables 
!
      real, dimension(:,:,:), allocatable:: DO
      real, dimension(:,:,:), allocatable:: BOD
      real, dimension(:,:,:), allocatable:: TEMP
      real, dimension(:,:,:), allocatable:: PO4
      real, dimension(:,:,:), allocatable:: P_Org
      real, dimension(:,:,:), allocatable:: NO2
      real, dimension(:,:,:), allocatable:: NO3
      real, dimension(:,:,:), allocatable:: NH4
      real, dimension(:,:,:), allocatable:: pH
      real, dimension(:,:,:), allocatable:: H2CO3
      real, dimension(:,:,:), allocatable:: HCO3
      real, dimension(:,:,:), allocatable:: CO3
      real, dimension(:,:,:), allocatable:: ALK
      real, dimension(:,:,:), allocatable:: CHLRE
      real, dimension(:,:,:), allocatable:: ALGAE_1
      real, dimension(:,:,:), allocatable:: ALGAE_2
      real, dimension(:,:,:), allocatable:: ZOO_1
      real, dimension(:,:,:), allocatable:: ZOO_2
      real, dimension(:,:,:), allocatable:: TDS
!
! Dimensioned and allocated water quality loading variables 
!
      real, dimension(:,:,:), allocatable:: TEMP_LOAD
      real, dimension(:), allocatable:: TEMP_TRIB_LOAD
!
      real, dimension(:,:,:), allocatable:: DO_LOAD
      real, dimension(:), allocatable:: DO_TRIB_LOAD
!
      real, dimension(:,:,:), allocatable:: BOD_LOAD
      real, dimension(:), allocatable:: BOD_TRIB_LOAD
!
      real, dimension(:,:,:), allocatable:: TDS_LOAD
      real, dimension(:), allocatable:: TDS_TRIB_LOAD
      
      
! 
! Tributary input
!
      real, dimension(:), allocatable:: DO_trib
      real, dimension(:), allocatable:: BOD_trib
      real, dimension(:), allocatable:: TEMP_trib
      real, dimension(:), allocatable:: PO4_trib
      real, dimension(:), allocatable:: P_Org_trib
      real, dimension(:), allocatable:: NO2_trib
      real, dimension(:), allocatable:: NO3_trib
      real, dimension(:), allocatable:: NH4_trib
      real, dimension(:), allocatable:: H2CO3_trib
      real, dimension(:), allocatable:: HCO3_trib
      real, dimension(:), allocatable:: CO3_trib
      real, dimension(:), allocatable:: ALK_trib
      real, dimension(:), allocatable:: CHLR_trib
      real, dimension(:), allocatable:: ALGAE_1_trib
      real, dimension(:), allocatable:: ALGAE_2_trib
      real, dimension(:), allocatable:: ZOO_1_trib
      real, dimension(:), allocatable:: ZOO_2_trib
      real, dimension(:), allocatable:: TDS_trib
!
! Nonpoint source concentrations
!
      real, dimension(:,:), allocatable:: DO_nps
      real, dimension(:,:), allocatable:: BOD_nps
      real, dimension(:,:), allocatable:: TEMP_nps
      real, dimension(:,:), allocatable:: PO4_nps
      real, dimension(:,:), allocatable:: P_Org_nps
      real, dimension(:,:), allocatable:: NO2_nps
      real, dimension(:,:), allocatable:: NO3_nps
      real, dimension(:,:), allocatable:: NH4_nps
      real, dimension(:,:), allocatable:: H2CO3_nps
      real, dimension(:,:), allocatable:: HCO3_nps
      real, dimension(:,:), allocatable:: CO3_nps
      real, dimension(:,:), allocatable:: ALK_nps
      real, dimension(:,:), allocatable:: CHLR_nps
      real, dimension(:,:), allocatable:: ALGAE_1_nps
      real, dimension(:,:), allocatable:: ALGAE_2_nps
      real, dimension(:,:), allocatable:: ZOO_1_nps
      real, dimension(:,:), allocatable:: ZOO_2_nps
      real, dimension(:,:), allocatable:: TDS_npsemp
!
! Point sources
!
      real, dimension(:,:), allocatable:: DO_source
      real, dimension(:,:), allocatable:: BOD_source
      real, dimension(:,:), allocatable:: TEMP_source
      real, dimension(:,:), allocatable:: PO4_source
      real, dimension(:,:), allocatable:: P_Org_source
      real, dimension(:,:), allocatable:: NO2_source
      real, dimension(:,:), allocatable:: NO3_source
      real, dimension(:,:), allocatable:: NH4_source
      real, dimension(:,:), allocatable:: TDS_source

end module block_wq
