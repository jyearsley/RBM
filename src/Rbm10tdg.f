c $Header: /home/CRTASS/columbia/temp_model/RCS/rbm10.f,v 2.1 2000/06/28 16:43:58 root Exp root $
c
c      PROGRAM RMAIN                                  
C
C     Dynamic river basin model for simulating water temperature
c     and total dissolved gas in the Columbia River system. 
c     This version uses Reverse Particle Tracking in the Lagrangian
c     mode and Lagrangian interpolation in the Eulerian mode
C     The Kalman filter is also implemented in this version
C     For additional information contact:
C
C     John Yearsley
C     EPA Region 10    ES-098
C     1200 Sixth Ave
C     Seattle, WA      98101
C     (206) 553-1532
C
      character*30 NAMEI
      INCLUDE 'rbm10tdg.fi'
C
C     Open file containing reach data
C
      WRITE(*,2600)
      read(*,1500) namei
      OPEN(UNIT=90,FILE=namei,STATUS='OLD')
c
c     Read header information from control file
c
      read(90,1500) namei
      write(*,2900) 
      write(*,1510) namei
c
c     Open file with hydrologic data
c
      open(unit=35,FILE=namei,STATUS='old')
C
C     Call systems programs to get started
C
C     SUBROUTINE BEGIN reads control file, sets up topology and 
C     important properties of reaches
C
      CALL BEGIN
C
C     SUBROUTINE SYSTMM performs the simulations
C
      CALL SYSTMM
C
C     Close files after simulation is complete
C
      CLOSE(UNIT=30)
      CLOSE(UNIT=35)
      CLOSE(UNIT=90)
 1500 FORMAT(A30)
 1510 format(1x,a30)
 1600 FORMAT(8F10.0)
 2600 FORMAT(' NAME OF FILE CONTAINING RIVER REACH DATA')
 2700 FORMAT(' NAME OF OUTPUT DATA FILE')
 2800 format(' Name of file with geometric data')
 2900 format(' Name of file with hydrologic data')
 3000 format(' Name of file with water quality data')
      STOP
      END
      SUBROUTINE BEGIN
      CHARACTER*5 type
      character*10 end_mark
      character*20 seg_name(200)
      character*30 NAMEI,plot_file
      integer start_date,end_date
	logical param(5)
      real*4 lat,long
      INCLUDE 'rbm10tdg.fi'
	equivalence (param(1),Temp)
      data nplot/200*0/,rm_trib/50*-99./,rm_plot/100*-99./
C
C     Initialize arrays and constants.  Filter parameters are set
C     for prediction mode (R_var(n)=0.0).  System variance (Q_var(n))
C     are set to values the same as estimates reported in the 
C     December 1999 Final Draft report of the Columbia River Temperature
C     Assessment
C
      Q_var(1,1)=0.007
      Q_var(2,1)=0.007
      Q_var(3,1)=0.002
      Q_var(4,1)=0.000
      Q_var(5,1)=0.000
      do n=1,5
	   do ncn=1,5
            R_var(n,ncn)=0.0
	   end do
      end do
      do n=1,5
         no_plots(n)=0
         do nn=-2,600
            trib(n,nn)=.FALSE.
	      do ncn=1,5
               P_var(n,nn,ncn,1)=Q_var(n,ncn)
               P_var(n,nn,ncn,2)=Q_var(n,ncn)
	      end do
         end do
         do nn=1,600
            type_res(n,nn)=.FALSE.
            type_riv(n,nn)=.FALSE.
	      do ncn=1,5
               Conc(n,nn,ncn,1)=0.0
	      end do
         end do
      end do
      no_inflow=0
      no_intrib=0
c
c     Total number of projects with spill is NO_PRJCT
c
	no_prjct=0
      nwprov=0
C
C     Card Group I
C               
c
c     Add for TDG
c
c
c     Number of constituents to be modeled (NO_CONST)
c
      no_const=0
      read(90,1150) param
 1150 format(20l1)
	do nprm=1,10
	   if(param(nprm)) no_const=no_const
	end do
c
      read(90,1200) start_date,end_date
      nyear1=start_date/10000
      nyear2=end_date/10000
      nysim=nyear2-nyear1+1
      ysim=nysim
      write(*,1200) start_date,end_date
 1200 format(/8i10)
      read(90,1200) no_rch
      write(*,1200) no_rch
c
 1400    format(16i5)
 1450    format(16(2x,a3))
C                                        
C     Card Group IIb. Reach characteristics
C
      do nr=1,no_rch
         write(*,*) ' Starting to read reach ',nr
         read(90,1480) rch_name(nr),no_elm(nr)
     .        ,no_plots(nr)
C
C     Setup locations for output.  Locations are specified in each
C     Reach by River Mile now, rather than by element number (NC)
C
         if(no_plots(nr).gt.0) then
            nplts=no_plots(nr)
            read(90,1065) (rm_plot(nr,np),np=1,nplts)
         end if
 1480    format(a20,6i10)
c
c     NC is a temporary element counter
c
         nc=0
c
c     NPRJ is the reach project counter
c
	   nprj=0
c
c     JP is a temporary plot location counter
c
         jp=1
c
c     Reading Reach Element information
c
         do ne=1,no_elm(nr)
            write(*,*) no_elm(nr),ne
C     
C     Card Type 3. Reach description, Reach type, begin and 
C     end river mile
C  
            read(90,1600) seg_name(ne),type,rmile1(nr,ne),rmile2(nr,ne)
            write(*,1600) seg_name(ne),type,rmile1(nr,ne),rmile2(nr,ne)
 1600       format(a20,5x,a5,10x,3f5.0)
c     
c     Read number of computational elements, weather province, 
c     headwaters number, number of entering tributaries, Reach number
c     if the tributary is one for which temperatures are simulated
c     
            read(90,1400) no_xsctn,kwtype,khead
     .           ,ntribs,nr_trib,prj_ndx
            write(*,1400) no_xsctn,kwtype
            x_dist1(nr,nc+1)=5280.*rmile1(nr,ne)
            nc1=nc+1
c
c     Reservoir segments
c
            if(type.eq.'RSRVR') then
               read(90,1700)
     .              rel_vol,rel_area,del_x
               write(*,1700) 
     .              rel_vol,rel_area,del_x
               do nx=1,no_xsctn
                  nc=nc+1
                  x_sctn=no_xsctn
                  type_res(nr,nc)=.TRUE.
c
c     Establish weather province
c
                  nwtype(nr,nc)=kwtype
                  if(kwtype.gt.nwprov) nwprov=kwtype
C     
C     Reservoir reaches
c     
                  dx(nr,nc)=5280.*del_x/x_sctn
c
c     Cross-sectional area and volume
c
                  res_area(nr,nc)=43560.*rel_area/x_sctn
                  res_vol(nr,nc)=43560.*rel_vol/x_sctn
c
c     Endpoints for reservoir segments
c
                  x_dist2(nr,nc)=x_dist1(nr,nc)-dx(nr,nc)
                  rmp_ft=5280.*rm_plot(nr,jp)
                  if(rmp_ft.ge.x_dist2(nr,nc).and.   
     .                 rmp_ft.lt.x_dist1(nr,nc)) then
                     nc_plot(nr,jp)=nc
                     jp=jp+1
                  end if
                  x_dist1(nr,nc+1)=x_dist2(nr,nc)
 1700             format(8f10.0)
               end do
c
c     Add for TDG
c
	         if(prj_ndx.gt.0) then
	            nprj=nprj+1
	            no_prjct=no_prjct+1
	            TDG_input(nr,nprj)=nc
                  TDG_model(nr,nprj)=prj_ndx
               end if
c			    
            end if
c     
c     River reaches
c     
            if(type.eq.'RIVER') then
               if(kwtype.gt.nwprov) nwprov=kwtype
               xsctn=no_xsctn
               delta_x=5280.*(rmile1(nr,ne)-rmile2(nr,ne))
               delta_x=delta_x/xsctn
               x_dist1(nr,nc+1)=5280.*rmile1(nr,ne)
               read(90,1700) a_a,b_a,a_w,b_w
               do nx=1,no_xsctn                
                  nc=nc+1
c
c     Establish weather province
c
                  nwtype(nr,nc)=kwtype
                  type_riv(nr,nc)=.TRUE.
                  a_area(nr,nc)=a_a
                  b_area(nr,nc)=b_a
                  a_width(nr,nc)=a_w
                  b_width(nr,nc)=b_w
                  x_dist2(nr,nc)=x_dist1(nr,nc)-delta_x
                  rmp_ft=5280.*rm_plot(nr,jp)
                  if(rmp_ft.ge.x_dist2(nr,nc).and.   
     .                 rmp_ft.lt.x_dist1(nr,nc)) then
                     nc_plot(nr,jp)=nc
                     write(*,*) 'river',nr,nc,jp,nc_plot(nr,jp)
                     jp=jp+1
                  end if
                  x_dist1(nr,nc+1)=x_dist2(nr,nc)
                  dx(nr,nc)=delta_x
               end do
            end if
            nc2=nc
c
c     Check to see if there is a tributary in this computational element.
c     If so, define pointers for later use.
c
            if(ntribs.gt.0) then
               no_inflow=no_inflow+1
               if(nr_trib.eq.0) then 
                  no_intrib=no_intrib+1
               else
                  trib_id(nr_trib)=no_inflow
               end if
               read(90,*)
               read(90,1850) trib_jnc
     .              ,rm_trib(no_inflow)
 1850          format(15x,i5,5x,f10.0)
               x_trib(no_inflow)=5280.*rm_trib(no_inflow)
               x_t=5280.*rm_trib(no_inflow)+1.
               do nnc=nc1,nc2
                  xdist1=x_dist1(nr,nnc)
                  xdist2=x_dist2(nr,nnc)
                  if(x_t.le.x_dist1(nr,nnc)
     .                 .and.x_t.gt.x_dist2(nr,nnc)) then
                     trib(nr,nnc)=.TRUE.
                     trib_ndx(nr,nnc)=no_inflow
                  end if
               end do
            end if
            read(90,1851) end_mark
 1851       format(a10)
            write(*,2851) end_mark
c            if(end_mark(1:3).ne.'End') pause 
 2851       format(' end_mark: ',a10)
C     
         end do
         x_dist2(nr,0) = x_dist1(nr,1)
         x_dist2(nr,-1)= x_dist2(nr,0) +dx(nr,1)
         x_dist2(nr,-2)= x_dist2(nr,-1)+dx(nr,1)
         x_dist2(nr,nc+1)=x_dist2(nr,nc)-dx(nr,nc)
         no_celm(nr)=nc
         do nnc=1,no_celm(nr)
            rm1=x_dist1(nr,nnc)/5280.
            rm2=x_dist2(nr,nnc)/5280.
         end do
      end do    
c 
c     Open meteorological file
c
      write(*,*) ' nwprov = ',nwprov
      DO  nw=1,NWPROV
         NWTAPE=10+nw
         read(90,1030) namei,evrate(nw)
         write(*,*) 'evrate = ',evrate(nw),' nw = ',nw
c         pause
         write(*,2700) 
c     
c     Open file with heat budget information
c     
         OPEN(UNIT=NWTAPE,FILE=NAMEI,STATUS='OLD')
         READ(NWTAPE,1025) WPNAME(nw)
         READ(NWTAPE,1027) nwpd,selev,lat,long,nwstrt,nwstop
c
c     Computational time step
c
         nw_year1=nwstrt/10000
         write(*,*) ' nwyear = ',nw_year1,nyear1,dt_comp
      end do
      xwpd=nwpd
      dt_comp=86400./xwpd
c     
c     Read header on advection file and advance to starting point
c     if necessary
c     
      read(35,1042) ny_adv1,ny_adv2
      ny_adv1=ny_adv1/10000
      if(ny_adv1.lt.nyear1) then
         do ny_adv=1,nyear1-ny_adv1
            do nd_adv=1,365
               read(35,*)
               if(tdg) then
                  do npj=1,no_prjct
                     read(35,*)
                  end do
	         end if
               do nt_adv=1,no_intrib
                  read(35,*)
               end do
            end do
         end do
      end if
c     
c     Set DT of boundary element equal to computational DT
c
      do nr=1,no_rch
         dt(nr,0)=dt_comp
      end do
 1020 FORMAT(A80)
 1025 format(a50)
 1027 format(5x,i5,3(5x,f5.0),2i10)
 1028 format(i5,7f10.0)
 1030 format(a30,f10.0)
 1035 format(1x,a30,e15.5)
 1040 FORMAT(8F10.0)
 1042 FORMAT(8I10)
 1043 format(i5,2f10.0,3i5)
 1044 FORMAT(16I5)
 1045 format(i5,f10.0,a60)
 1048 FORMAT(8F10.0)
 1050 FORMAT((A30,5X,A5,10x,5F5.0))
 1060 FORMAT(2F10.0,A20/16F5.0)
 1063 FORMAT(F10.0,A20/16F5.0)
 1065 FORMAT(16F5.0)
 1080 FORMAT(A3)
 1085 format(1x,a3)
 1145 FORMAT(8F10.2)
 1152 FORMAT(6I3)
 1250 format(2i5,3f10.0)
 1500 format(i9)
 2500 FORMAT(' ENERGY BUDGET FILE FOR METEOROLOGIC PROVINCE - ',I5)
 2700 format(' energy budget file')
 3000 FORMAT(1H0,'CARD SEQUENCE ERROR IN DATA FOR REACH - ',I5)
 3500 format(' reservoir flow file')
C 
C     ******************************************************
C                         Return to RMAIN
C     ******************************************************
C
      RETURN
  900 END
      SUBROUTINE SYSTMM
      real*4 WDATA(5,7),EDATA(7),xa(4),ta(4),var(4)
     .      ,dt_part(600),x_part(600),C0(5),Var0(5)
     .      ,q_spill(20)
      integer no_dt(600),simyr,strt_elem(600)
     .     ,ndltp(4),nterp(4)
      real*4 plot_data(200,2)
      INCLUDE 'rbm10tdg.fi'
      EQUIVALENCE (EDATA(1),QNS)
      data ndltp/-2,-1,-2,-2/,nterp/4,3,2,3/
      data pi/3.14159/,rfac/304.8/
C
      Diff_TDG=2.0e-05/(12*12*2.54*2.54)
      time=0.0
      n1=1
      n2=2
      nobs=0
      simyr=0
 50   continue
      simyr=simyr+1
      write(*,*) ' Simulation Year - ',simyr
      DO ND=1,365
         nobs=nobs+1
         DO NDD=1,NWPD
C     
C     Read weather data from files if time period is correct
C     
            do nw=1,nwprov
               nwr=10+nw
               READ(NWR,1028) LDUMM,(WDATA(nw,nnw),nnw=1,7)
            end do
c     
c     begin reach computations
c     
            ind=nd
            ipd=ndd    
            day=nd
c     
            if(ndd.ne.1) go to 90
c     
c     Read flow and water quality
c     
c     Main stem inflows and outflows for each reach first
c     Flows are cumulative and do not include tributaries if
c     tributaries are downstream of the inflow junction
c     
c     
c     Headwaters flow and temperature
c   
            pi=3.14159
            ttme=nd
	      do nr=1,no_rch
               read(35,1550) jy,jobs
     .                  ,(qin(nr,1),C_head(nr,ncn),ncn=1,no_rch)
	      end do
            if(TDG) read(35,1600) (nppj,q_spill(npj),npj=1,no_prjct)
        
c     
c     Check for tributaries
c     
            if(no_intrib.gt.0) then
c     
c     Tributary flow and temperature (input)
c     
               do in=1,no_intrib
                  read(35,1600)
     .                 in_trb,q_trib(in_trb)
     .                ,(C_trib(in_trb,ncn),ncn=1,no_const)
c
c     Remove comment ("c" in Column 1) for Scenario 3
c
c                  if(T_trib(in_trb).gt.16.0) T_trib(in_trb)=16.0
               end do
            end if
c     
c     Tributary temperatures (computed)
c     
            if(no_rch.gt.1) then
               do nr=1,no_rch-1
                  nt=trib_id(nr)
                  do ncn=1,no_const
	               C_trib(nt,1)=Conc(nr,no_celm(nr),ncn,n1)
				end do 
               end do
            end if
c
c     Call the flow balance subroutine to set up the
c     system hydraulics
c
            call balance
c
 1400       format(4i5,4f10.0)
 1500       format(i5,i10,10F10.0)
 1510       format(8f10.5)
 1550       format(2i5,10f10.0)
 1600       format(i5,f10.0,f5.0)
 90         continue
            nsmpl=1
c     
c     Begin cycling through the reaches
c     
            do nr=1,no_rch
               nc=no_celm(nr)
               qsum=qin(nr,1)
	         do ndmm=-2,0
	            do ncn=1,no_const
                     conc(nr,ndmm,ncn,n1)=C_head(nr,1)
	            end do
               end do
               do ncn=1,no_const
	            conc(nr,nc+1,ncn,n1)=C_head(nr,ncn)
                  P_var(nr,nc+1,ncn,n1)=P_var(nr,nc,ncn,n1)
	         end do
               x_head=x_dist1(nr,1)
               x_bndry=x_head-5.0
c     
c     First do the reverse particle tracking
c     
               do nc=no_celm(nr),1,-1
                  nx_s=1
                  nx_part=nc
c                  dt_part(nc,nx_s)=dt(nr,nx_part)
c                  dt_total=dt_part(nc,nx_s)
                  dt_part(nc)=dt(nr,nx_part)
                  dt_total=dt_part(nc)
                  x_part(nc)=x_dist2(nr,nx_part)
 100              continue
c
c     Determine if the total elapsed travel time is equal to the
c     computational interval
c
                  if(dt_total.lt.dt_comp) then
                     x_part(nc)=x_part(nc)+dx(nr,nx_part)
c
c     If the particle has started upstream from the boundary point, give it
c     the value of the boundary
c
                     if(x_part(nc).ge.x_bndry) then
                        x_part(nc)=x_head
                        dt_part(nc)=dt(nr,1)
                        go to 200
                     end if
c
c     Increment the segment counter if the total time is less than the
c     computational interval
c
                     nx_s=nx_s+1
                     nx_part=nx_part-1
                     dt_part(nc)=dt(nr,nx_part)
                     dt_total=dt_total+dt_part(nc)
                     go to 100
                  else
c
c     For the last segment of particle travel, adjust the particle location
c     such that the total particle travel time is equal to the computational 
c     interval.
c
                     dt_part(nc)
     .                    =dt_comp-dt_total+dt_part(nc)
                     x_part(nc)=x_part(nc)
     .                    +u(nr,nx_part)*dt_part(nc)
                     if(x_part(nc).ge.x_head) then
                        x_part(nc)=x_head
                        nx_s=nx_s-1
                        dt_part(nc)=dt(nr,1)
                     end if
 3700                format(f10.1)
 3710                format(2i5,3f10.1)
                  end if
 200              continue
                  if(nx_part.lt.1) nx_part=1
                  strt_elem(nc)=nx_part
                  no_dt(nc)=nx_s
               end do
               do nc=1,no_celm(nr)
c     
c     read meteorological data from the appropriate file
c     
                  iwr=nwtype(nr,nc)
                  if(iwr.eq.0) go to 250
c    
c     Correspondence table for energy budget terms
c     
c     Net solar radiation (kcal/meter^2/second)
c     
c     qns=wdata(iwr,1)
c     
c     Net atmospheric radiation (kcal/meter^2/second)
c     
c     qna=wdata(iwr,2)
c     
c     Dry bulb temperature (deg C)
c     
c     dbt=wdata(iwr,3)
c     
c     Wind speed (meters/second)
c     
c     wind=wdata(iwr,4)
c     
c     Factor for Bowen ratio ((deg C)^-1)
c     
c     pf=wdata(iwr,5) 
c     
c     Vapor pressure at given air temperature (mb)
c     
c     ea=wdata(iwr,6)
c     
c     Photo period (fraction of a day.  Not used in the energy budget)
c     
c     phper=wdata(iwr,7)
c     
                  do nnw=1,7
                     edata(nnw)=wdata(iwr,nnw)
                  end do
 250              continue
c     
c     Now do the third-order interpolation to
c     establish the starting temperature values
c     for each parcel
c     
                  ncell=strt_elem(nc)
                  npndx=1
                  do ntrp=2,4
                     ntest=ncell+ntrp-3
                     if(trib(nr,ntest)) then
                        npndx=ntrp
                     end if
                  end do
c                  ndltp=-2
c     
c     If starting element is the first one, then set
c     the initial temperature to the boundary value
c     
                  if(ncell.eq.1) then
	               do ncn=1,nconst
	                  C0(ncn)=C_head(nr,ncn)
	               end do
                     go to 350
                  end if
c     
c     Perform polynomial interpolation
c     
c
c     Add for TDG
c
                  do ncn=1,no_const                 
                     do ntrp=1,nterp(npndx)
                        npart=ncell+ntrp+ndltp(npndx)-1
                        xa(ntrp)=x_dist2(nr,npart)
                        ta(ntrp)=Conc(nr,npart,ncn,n1)
                        var(ntrp)=P_var(nr,npart,ncn,n1)
                     end do
                     x=x_part(nc)
c     
c     Call the interpolation function
c
                     C0(ncn)=tntrp(xa,ta,x,nterp(npndx))
                     Var0(ncn)=tntrp(xa,var,x,nterp(npndx))
	            end do
 300              continue
 350              continue
                  dt_calc=dt_part(nc)
                  do nm=no_dt(nc),1,-1
                     nw=nwtype(nr,ncell)
	               vel=u(nr,ncell)
                     z=depth(nr,ncell)
c
c     Transfer of heat across the air-water interface
c
                        call energy(C0(1),qsurf,A,B,nw)
                        qdot=qsurf/(z*rfac)
                        C0(1)=C0(1)+qdot*dt_calc
                        if(C0(1).lt.0.0) C0(1)=0.0
                        phi=(1.+A/(z*rfac)*dt_calc)
                        var0(1)=phi*var0(1)*phi+Q_var(nr,1)
c
c     Transfer of total dissolved gas across the air-water interface
c
	                  if(tdg) then
				         TDG_sat=21.1-0.3125*Conc(nr,ncell,1,n1)
	                     K_TDG=sqrt(Diff_TDG*vel/z)
	                     C0(2)=C0(2)+K_TDG*(TDG_sat-C0(2))*dt_calc
	                  end if
 400                 continue
c
c     Look for a tributary.  If none, skip to STATEMENT 450
c
                     if(.not.trib(nr,ncell)) go to 450
                     q1=qin(nr,ncell)
                     q2=qin(nr,ncell+1)
                     nt_trb=trib_ndx(nr,ncell)
	               do ncn=1,no_const
                        C0(ncn)=(q1*C0(ncn)
     .	                      +q_trib(nt_trb)*C_trib(nt_trb,ncn))
     .                          /q2
	               end do
  450                continue
c
c    If TDG is being simulated, check to see if this element
c    contains  project with TDG input
c
                     if(TDG) then
	               if(TDG_input(nr,nprj).eq.ncell) then
	                 TDGF=C0(2)
	                 QS=q_spill(nprj)
	                 QT=q2
	                 model=TDG_model(nr,nprj)
	                 go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
     .                        model
	                 go to 20   
    1 				    C0(2)=C0(2)*dwr_tdg(QS,QT,TDGF)
                       go to 20	   
    2					    C0(2)=C0(2)*lwg_tdg(QS,QT,TDGF)
                       go to 20
    3					    C0(2)=C0(2)*lgs_tdg(QS,QT,TDGF)
                       go to 20
    4					    C0(2)=C0(2)*lmn_tdg(QS,QT,TDGF)
                       go to 20
    5					    C0(2)=C0(2)*ihr_tdg(QS,QT,TDGF)
                       go to 20
    6					    C0(2)=C0(2)*chj_tdg(QS,QT,TDGF)
                       go to 20
    7					    C0(2)=C0(2)*wel_tdg(QS,QT,TDGF)
                       go to 20
    8					    C0(2)=C0(2)*ris_tdg(QS,QT,TDGF)
                       go to 20
    9					    C0(2)=C0(2)*rrh_tdg(QS,QT,TDGF)
                       go to 20
   10					    C0(2)=C0(2)*wan_tdg(QS,QT,TDGF)
                       go to 20
   11					    C0(2)=C0(2)*prd_tdg(QS,QT,TDGF)
                       go to 20
   12					    C0(2)=C0(2)*mcn_tdg(QS,QT,TDGF)
                       go to 20
   13					    C0(2)=C0(2)*jda_tdg(QS,QT,TDGF)
                       go to 20
   14					    C0(2)=C0(2)*tda_tdg(QS,QT,TDGF)
                       go to 20
   15					    C0(2)=C0(2)*bon_tdg(QS,QT,TDGF)
                       go to 20
   20                  continue
                       nprj=nprj+1
				   end if
	               end if
c  
c    Update constituent values
c
                     do ncn=1,no_const
                        Conc(nr,nc,ncn,n2)=C0(ncn)
                        P_var(nr,nc,ncn,n2)=Var0(ncn)
	               end do
                     ncell=ncell+1
                     dt_calc=dt(nr,ncell)
                  end do
               end do
               
 260           continue
c
c     Set up the output
c
               do np=1,no_plots(nr)
                  nc_plt=nc_plot(nr,np)
                  plot_data(np,1)=Conc(nr,nc_plt,1,n1)
                  plot_data(np,2)=Conc(nr,nc_plt,2,n1)
c                  plot_data(np,2)=0.0
               end do
               time=nd
               xdd=ndd
               clock=(xdd-0.5)*dt_comp
               time=time+clock/86400.
               time=nyear1+simyr-1.0+time/365.
               zero=0.0
               if(no_plots(nr).gt.0) then
                  itplot=40+nr
                  write(itplot,4700) time,(rm_plot(nr,np)
     .                 ,plot_data(np,1),plot_data(np,2)
     .                 ,np=1,no_plots(nr))
               end if
            end do
            ntmp=n1
            n1=n2
            n2=ntmp
c     
c     End of weather period loop
c     
 4700       format(f10.4,8(2f6.1,f6.3))
 4750       format(f10.4,10(i4,f8.0))
         end do
C     
c     End of main loop
c     
      end do
c     
c     Check if there are years remaining to be simulated.  If so,
c     start at the top (STATEMENT 50)
c
      if(simyr.lt.nysim) go to 50
c     
c     FORMAT statements
c     
 1020 FORMAT(A80)
 1025 format(a50)
 1027 format(5x,i5,3(5x,f5.0),2i10)
 1028 format(i5,7f10.0)
 1042 format(8i10)
 1045 format(i5,f10.0,a60)
c 
c     ******************************************************
c                        return to rmain
c     ******************************************************
c
  950 return
      end
      subroutine balance
      include 'rbm10tdg.fi'
      dt_max=-10000.
      do nr=1,no_rch
         do nc=1,no_celm(nr)
            qsum=qin(nr,nc)
            if(type_res(nr,nc)) then
               s_area=res_area(nr,nc)
               volume=res_vol(nr,nc)
               depth(nr,nc)=volume/s_area
               dt(nr,nc)=volume/qin(nr,nc)
               u(nr,nc)=dx(nr,nc)/dt(nr,nc)
            end if
            if(type_riv(nr,nc)) then
               x_area=a_area(nr,nc)*(qsum**b_area(nr,nc))
               x_wide=a_width(nr,nc)*(qsum**b_width(nr,nc))
               depth(nr,nc)=x_area/x_wide
               u(nr,nc)=qsum/x_area
               dt(nr,nc)=dx(nr,nc)/u(nr,nc)
            end if
            if(trib(nr,nc)) then
               nt_trb=trib_ndx(nr,nc)
               qsum=qsum+q_trib(nt_trb)
            end if
            if(dt(nr,nc).gt.dt_comp) then
               nmax=nc
               dt_max=dt(nr,nc)
            end if
            qin(nr,nc+1)=qsum
         end do
         dt(nr,no_celm(nr)+1)=dt_comp
         nt_trb=trib_id(nr)
         q_trib(nt_trb)=qsum
      end do
      return
      end
      SUBROUTINE ENERGY(TSURF,QSURF,A,B,NW)
      REAL*4 LVP
      real*4 q_fit(2),T_fit(2)
      INCLUDE 'rbm10tdg.fi'
      T_fit(1)=tsurf-1.0
      T_fit(2)=tsurf+1.0
      do i=1,2
         E0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
         RB=PF*(DBT-T_fit(i))
         LVP=597.0-0.57*T_fit(i)
         QEVAP=1000.*LVP*EVRATE(NW)*WIND
         if(qevap.lt.0.0) qevap=0.0
         QCONV=RB*QEVAP
         QEVAP=QEVAP*(E0-EA)
         QWS=6.693E-2+1.471E-3*T_fit(i)
         q_fit(i)=QNS+QNA-QWS-QEVAP+QCONV
      end do
c
c     q=AT+B
c
      A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
      B=(T_fit(1)*q_fit(2)-T_fit(2)*q_fit(1))
     .     /(T_fit(1)-T_fit(2))
      qsurf=0.5*(q_fit(1)+q_fit(2))
C 
C     ******************************************************
C               Return to Subroutine RIVMOD
C     ******************************************************
C
      RETURN
      END
      function nodays(jtime,jy0)
      dimension ndmo(12)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334/
      jy=(jtime/10000)
      jrem=(jtime-jy*10000)
      jm=jrem/100
      jd=jrem-jm*100 
      ny=jy-jy0
      nodays=365*ny+ndmo(jm)+jd
      return
      end 
      function ndate(n,ny0)
      dimension ndmo(13)
      data ndmo/0,31,59,90,120,151,181,212,243,273,304,334,365/
      ny=(n-1)/365
      njul=n-ny*365
      nm=1
 10   continue
      if(ndmo(nm+1).ge.njul) go to 50
      nm=nm+1
      go to 10
 50   continue
      nday=njul-ndmo(nm)
      nmon=100*nm
      nyear=10000*(ny0+ny)
      ndate=nyear+nmon+nday
      return
      end
c
c	Third-order polynomial interpolation using Lagrange
c     polynomials.  FUNCTION is SUBROUTINE POLINT from
c     Numerial Recipes
c
      FUNCTION tntrp(XA,YA,X,n)
c      PARAMETER (N=4) 
c      DIMENSION XA(N),YA(N),C(N),D(N)
      DIMENSION XA(4),YA(4),C(4),D(4)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) DEN=0.001
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
	tntrp=y
      RETURN
      END
      function dwr_tdg(qs,qt,tdgf)
      real*8 x1,x2,netsum
      x1=2.*qs/9.3-1.
      x2=2.*(qt-1.3)/17.8-1.
      netsum=-2.012780008392292
      netsum=netsum-3.675795149953003*x1
      netsum=netsum+9.524299666977761*x2
      netsum=netsum-1.924227033305504*x1*x1
      netsum=netsum-6.624543997979434*x2*x2
      netsum=netsum+0.4752812638142955*x1*x1*x1
      netsum=netsum-5.859428239806556*x2*x2
      netsum=netsum+11.07944943048697*x1*x2
      dwr_tdg=20.6*(netsum+1)/2.+98.1
      return
      end
      function lwg_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/84.6-1.
      x2=2.*(qt-17.1)/170.4-1.
      x3=2.*(tdgf-96.2)/13.8-1.
      netsum=0.1588596107535112
      netsum=netsum+0.5744958427905256*x1
      netsum=netsum+0.3578651466760778*x2
      netsum=netsum-0.2193197501447088*x1*x3
      lwg_tdg=30.4*(netsum+1)/2.+96.5
      return
      end
      function lmn_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/66.1-1.
      x2=2.*(qt-14.5)/171.6-1.
      x3=2.*(tdgf-99.1)/27.6-1.
      netsum=0.4592529537481663
      netsum=netsum+0.5379049287564071*x1
      netsum=netsum+0.1439123744998875*x2
      netsum=netsum+0.4576569580011143*x3
      netsum=netsum-0.3897640533009928*x2*x2
      lmn_tdg=33.3*(netsum+1)/2.+95.4
      return
      end
      function lgs_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/73.6-1.
      x2=2.*(qt-13.9)/164.6-1.
      x3=2.*(tdgf-95.8)/28.3-1.
      netsum=0.6878084468686018
      netsum=netsum+1.14279250239962*x1
      netsum=netsum-0.4746525782098229*x2
      netsum=netsum+0.4745354908363749*x3
      netsum=netsum+0.3283477952191652*x1*x1
      netsum=netsum-0.9416203212620867*x1*x2
      lgs_tdg=27.2*(netsum+1)/2.+96.2
      return
      end
      function ihr_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/125.9-1.
      x2=2.*(qt-13.)/182.9-1.
      x3=2.*(tdgf-95.6)/29.7-1.
      netsum=0.1341272677776274
      netsum=netsum+0.4485584619259884*x1
      netsum=netsum+0.2182546160494846*x2
      netsum=netsum-0.7826736889180375*x1*x3
      netsum=netsum+0.695763675289968*x2*x3
      netsum=netsum+0.3223848290352992*x1*x2*x3
      ihr_tdg=28.8*(netsum+1)/2.+98.3
      return
      end
      function chj_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/23.6-1.
      x2=2.*(qt-68.)/113.6-1.
      x3=2.*(tdgf-102.9)/12.6-1.
      netsum=0.3485352887924531
      netsum=netsum+0.528236037228155*x1
      netsum=netsum+0.6065521585482949*x3
      netsum=netsum-9.156290532428719D-002*x3*x3
      netsum=netsum-0.1340116712038874*x2*x3
      chj_tdg=17.4*(netsum+1)/2.+102.9
      return
      end
      function wel_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/54.3-1.
      x2=2.*(qt-67.6)/139.1-1.
      x3=2.*(tdgf-102.4)/11.3-1.
      netsum=-0.1594964567638874
      netsum=netsum+0.2307011514354815*x2
      netsum=netsum+0.7433482537124089*x3
      netsum=netsum+0.3400283946461143*x1*x1*x1
      netsum=netsum-0.2763414145641050*x3*x3*x3
      netsum=netsum+0.3361595273322946*x1*x2
      wel_tdg=16.8*(netsum+1)/2.+103.0
      return
      end
      function ris_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/169.0-1.
      x2=2.*(qt-45.8)/293.2-1.
      x3=2.*(tdgf-99.5)/30.7-1.
      netsum=0.5205350181937786
      netsum=netsum+0.6792921635306042*x1
      netsum=netsum-0.4765834853683446*x2
      netsum=netsum+0.5578479593803300*x3
      netsum=netsum-0.5319871634773443*x3*x3
      netsum=netsum+0.2070421416267055*x2*x2*x2
      netsum=netsum-0.7551993772999508*x1*x2
      netsum=netsum+0.7757959962032884*x2*x3
      dwr_tdg=32.60001*(netsum+1)/2.+102.7
      return
      end
      function rrh_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/200.6-1.
      x2=2.*(qt-63.5)/277.2-1.
      x3=2.*(tdgf-100.8)/31.3-1.
      netsum=0.2293710452132839
      netsum=netsum+0.2829624874267437*x1
      netsum=netsum+0.3893788606986859*x3
      netsum=netsum-0.6067024145801769*x1*x3
      netsum=netsum+0.5033555918168181*x2*x3
      netsum=netsum+0.2433066263205654*x1*x2*x3
      rrh_tdg=31.4*(netsum+1)/2.+102.1
      return
      end
      function wan_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/69.9-1.
      x2=2.*(qt-72.)/152.9-1.
      x3=2.*(tdgf-101.4)/14.5-1.
      netsum=0.4275705474781085
      netsum=netsum+0.7916471670144922*x1
      netsum=netsum-0.2901296584860031*x2
      netsum=netsum+0.6704809013244161*x3
      netsum=netsum-0.4745101531408344*x2*x2
      wan_tdg=18.5*(netsum+1)/2.+101.5
      return
      end
      function prd_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/168.9-1.
      x2=2.*(qt-77.1)/151.6-1.
      x3=2.*(tdgf-100.4)/20.7-1.
      netsum=0.4857522925162594
      netsum=netsum+0.5887877594348601*x1
      netsum=netsum+0.4422592773477501*x3
      netsum=netsum-0.2224974137745099*x1*x1
      netsum=netsum-0.313001866590213*x1*x3
      dwr_tdg=20.4*(netsum+1)/2.+101.8
      return
      end
      function mcn_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/218.9-1.
      x2=2.*(qt-97.1)/291.5-1.
      x3=2.*(tdgf-99.7)/19.6-1.
      netsum=0.3780134616741254
      netsum=netsum+0.664990671302595*x1
      netsum=netsum+6.881348794175454D-002*x3
      netsum=netsum-0.1819254166756301*x1*x2
      netsum=netsum-0.1786138846673968*x1*x3
      netsum=netsum+0.3605859359333515*x1*x2*x3
      mcn_tdg=27.1*(netsum+1)/2.+99.6
      return
      end
      function jda_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/119.6-1.
      x2=2.*(qt-93.9)/297.3-1.
      x3=2.*(tdgf-100.6)/18.8-1.
      netsum=0.2937524242470778
      netsum=netsum+0.3379218902737073*x1
      netsum=netsum+0.3747977305210063*x2
      netsum=netsum+0.4724854650873372*x3
      netsum=netsum-0.3773630356169043*x1*x3
      netsum=netsum-0.2829591805303990*x1*x2*x3
      jda_tdg=25.1*(netsum+1)/2.+98.8
      return
      end
      function tda_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/407.3-1.
      x2=2.*(qt-87.2)/483.-1.
      x3=2.*(tdgf-82.4)/50.0-1.
      netsum=0.1829053614213454
      netsum=netsum+0.2786787677653985*x1
      netsum=netsum-0.3769435603276781*x2
      netsum=netsum+0.7394767667623935*x3
      netsum=netsum-0.3302600036176216*x1*x1
      netsum=netsum-0.3523829164363117*x2*x2
      netsum=netsum-0.5312171226484144*x2*x2*x2
      netsum=netsum-0.9317877150165886*x3*x3*x3
      netsum=netsum-0.4465446527441894*x1*x2
      netsum=netsum-0.6443473945850536*x1*x3
      netsum=netsum+1.6660574975925700*x2*x3
      netsum=netsum+1.8890572085047500*x1*x2*x3
      dwr_tdg=20.6*(netsum+1)/2.+98.1
      return
      end
      function bon_tdg(qs,qt,tdgf)
      real*8 x1,x2,x3,netsum
      x1=2.*qs/210.-1.
      x2=2.*(qt-98.4)/280.9-1.
	x3=2.*(tdgf-82.4)/33.1-1.
      netsum=5.145860526037795D-002
      netsum=netsum+0.8106134425929802*x1
      netsum=netsum+0.2989893389868148*x2
      netsum=netsum+1.1194864704397040*x3*x3
      netsum=netsum-0.6543762829781076*x1*x3
      netsum=netsum-0.7025646372586098*x2*x3
      bon_tdg=22.2*(netsum+1)/2.+98.8
      return
      end





