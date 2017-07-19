SUBROUTINE Water_Balance
!
USE Block_Hydro
USE Block_Network
!
IMPLICIT NONE
!
integer     :: nc,ncnf,nr,nrc,nrc_1,nrc_2,ntrib,nntrib,total_cells
integer     :: tot_trib_in,trib_seg_in
real        :: Q_start,Q_sum_trib
!
logical     :: TRIBS_DONE
!
total_cells=0
!
do nr = 1,nreach
! 
  do nc = 1,no_cells(nr)
    total_cells=total_cells+1
    nrc = reach_cell(nr,nc)
    if (nc .eq. 1) then
!
!  Flow into the headwaters cell is the same as the flow out
!
      Q_in(nrc) = Q_out(nrc)
    else
!
!  Flow into downstream cells is equal to the flow out of the most upstream cell
!
      Q_in(nrc) = Q_out(nrc-1)
    end if
    ntrib=no_tribs(nrc)
    Q_sum_trib = 0.0
    if (ntrib .gt. 0) then
      do nntrib  = 1,ntrib
        ncnf = conflnce(nrc,nntrib)
        Q_sum_trib = Q_sum_trib + Q_out(ncnf)
      end do
    end if
!
! Distributed flow is divided equally among each the segments in a cell
!
    Q_diff(nrc) = (Q_out(nrc) - Q_in(nrc) - Q_sum_trib)/ndelta(nrc)
!
!    
  end do
!
! Now define the inflow to each computational segment
!
!
  TRIBS_DONE = .FALSE.
!
! Initialize the inflow with the VIC inflow and cell number at the headwaters segment
!
  Q_in_seg(nr,1) = Q_in(reach_cell(nr,1))
  nrc_1 = segment_cell(nr,1)
!
  do nc = 1,no_celm(nr)-1
    nrc_2 = segment_cell(nr,nc)
!
! If this is a new cell, reset the TRIBS_DONE flag
!
    if (nrc_1 .ne. nrc_2) TRIBS_DONE = .FALSE.
    ntrib = no_tribs(nrc_2)
    Q_sum_trib = 0.0
    if (ntrib .gt. 0 .and. .not.TRIBS_DONE) then
    do nntrib = 1,ntrib
        ncnf = conflnce(nrc_2,nntrib)
        Q_sum_trib = Q_sum_trib + Q_out(ncnf)
    end do
    end if
    TRIBS_DONE = .TRUE.
    nrc_1 = nrc_2
    Q_out_seg(nr,nc) = Q_in_seg(nr,nc) + Q_sum_trib + Q_diff(nrc_2) 
    Q_in_seg(nr,nc+1) = Q_out_seg(nr,nc) 
  end do
!
! Outflow from the last segment (no_celm(nr)). The distributed flow is not 
! included in this since we're assuming here that the last segment is 
! really part of the confluent stream. (Food for thought)
!
  Q_out_seg(nr,no_celm(nr)) = Q_in_seg(nr,no_celm(nr))  
end do
end SUBROUTINE Water_Balance
