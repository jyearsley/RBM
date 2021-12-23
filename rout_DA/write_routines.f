c SUBROUTINES RELATED TO WRITING
c write_data()
c write_month()
c

      SUBROUTINE WRITE_DATA
c
c  Add Qin - JRY - 9/30/2009
c
c     & (FLOW, Qin, DAYS, NAME5, FACTOR_SUM, OUTPATH,IDAY,IMONTH,IYEAR)
     & (write_flow,rbm10,FLOW, Qin, DAYS, NAME5
     & ,FLOWOUT,fdleno,HEATPATH,hdleni,loc,clen
     & ,Flow_Cells,Force_Cells,Nflow,Nflow_pd,Nheat,Nheat_pd 
     & ,nday,IDAY,IMONTH,IYEAR,skip_MET
     & ,aa_d,bb_d,aa_w,bb_w)

      IMPLICIT NONE

      INTEGER DAYS
c
c  Add Qin - JRY - 9/30/2009
c
      REAL    FLOW(DAYS),Qin(days),Dmmy(7),Heat_data(7)
      INTEGER IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER I, ii, CLEN, FLEN, Flow_Cells, Force_Cells
     &       ,nday,Nflow,Nflow_pd, navg, Nheat, Nheat_pd, nrec
     &       ,fdleni,fdleno,hdleni,hdleno
     &       ,skip_MET,nyy,nmm,ndd,nhh
     &       ,Ntest
      integer dmmy_day
      logical rbm10,write_flow,write_header
      CHARACTER*24 NAME5   !was 5 6/25/2009
      REAL    FACTOR_SUM
     &    ,flowavg,flowin,Q_dmmy,depth,width,vel,heat_pd
     &    ,yy,mm,dd,hh
      real    aa_d,bb_d,aa_w,bb_w
c
      character*25 loc
      character*50 header
      CHARACTER*72 OUTPATH
      character*72 flowout,flowpath,heatout,heatpath
      write (*,*) 'NDAY - ',nday
      if (rbm10) then
         heat_pd=Nheat_pd
c
c     Open forcing file from VIC
c 
         open(25,FILE=heatpath(1:hdleni)//loc(1:clen),status='old')
c
c     Read header on forcing file
c 
        do ii=1,6
            read(25,'(a)') header
        end do
c
c     Skip records if necessary
c 
        if (skip_MET.gt.0) then
           do ii=1,skip_MET
              read(25,*)
           end do
         end if  
         do i=1,nday
           if (write_flow) then
              nrec=Flow_Cells*(i-1)+Nflow
              flowavg=0.5*(flowin+flow(i))
              depth=aa_d*(flowavg**bb_d)
              width=aa_w*(flowavg**bb_w)
              vel=flowavg/(depth*width)
c      
c     Write the flow for this grid cell to the RBM direct access file
c
             write(15,'(2i5,3f10.1,2f6.2)',rec=nrec)
     &              i,Nheat,Q_dmmy,flow(i),Q_dmmy,depth,vel
c
           end if
             do navg=1,7
                Heat_data(navg)=0.0
             end do
             do navg=1,Nheat_pd
c
c            Read full_data output from VIC runs and average
c            year,month,day,air_temp,vp,short_wave,long_wave
c     &                ,rho,press,wind
c
c               read(25,*) dmmy_day,dmmy
               read(25,*) nyy,nmm,ndd,nhh,dmmy
c Air temperature - deg C
               Heat_data(1)=Heat_data(1)+dmmy(1)/heat_pd
c Vapor pressure - kPa
               Heat_data(2)=Heat_data(2)+dmmy(2)/heat_pd
c Shortwave radiation Watts/m**2
               Heat_data(3)=Heat_data(3)+dmmy(3)/heat_pd
c Longwave radiation Watts/m**2
               Heat_data(4)=Heat_data(4)+dmmy(4)/heat_pd
c Air density Kg/m**3
               Heat_data(5)=Heat_data(5)+dmmy(5)/heat_pd
c Atmospheric pressure - kPa               
               Heat_data(6)=Heat_data(6)+dmmy(6)/heat_pd
c Wind speed - m/sec
               Heat_data(7)=Heat_data(7)+dmmy(7)/heat_pd
            end do
            nrec=Force_Cells*(i-1)+Nheat
c
c         Write the heat budget data for this grid cell
c         to the RB direct access file
c         
            write(*,*) nday,Heat_data
            write(16,'(2i5,f6.1,f6.3,2f7.1,f6.3,f7.1,f5.1)',rec=nrec)
     &           i,Nheat,Heat_data
     
         end do
         close(25)
      else 
         OPEN(30, FILE = flowout(1:fdleno)//loc(1:CLEN)//'.day'
     .    ,status='unknown')
         DO I = 1,DAYS
            WRITE(30,*) IYEAR(I),IMONTH(I),IDAY(I),FLOW(I),Qin(i)
         END DO
      end if
      CLOSE(30)
      RETURN
      END
