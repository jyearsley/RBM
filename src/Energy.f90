<<<<<<< HEAD
SUBROUTINE Energy(T_surf,q_surf,ncell,nr)
   use Block_Energy
   implicit none
   integer::i,ncell,nd,nr
   real::A,B,e0,q_surf,q_conv,q_evap,q_ws,td,T_surf
   real, dimension(2):: q_fit, T_fit
!
   td=nd
   T_fit(1)=T_surf-1.0
   T_fit(2)=T_surf+1.0
   do i=1,2
      e0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
      rb=pf*(dbt(ncell)-T_fit(i))
      lvp=597.0-0.57*T_fit(i)
      q_evap=1000.*lvp*evap_coeff*wind(ncell)
      if(q_evap.lt.0.0) q_evap=0.0
      q_conv=rb*q_evap
      q_evap=q_evap*(e0-ea(ncell))
      q_ws=6.693E-2+1.471E-3*T_fit(i)
      q_fit(i)=q_ns(ncell)+q_na(ncell)-q_ws-q_evap+q_conv
!if (nr .eq. 1) write(26,*) 'energy',ncell,T_surf,q_ws
   end do
!
!     q=AT+B
=======
SUBROUTINE Energy(T_surf,q_surf,ncell,nrr_tmp)
   use Block_Energy
   implicit none
   integer::i,ncell,nd,nrr_tmp
   real::A,B,e0,q_surf,q_conv,q_evap0,q_evap,q_ws,td,T_surf
   real::dbt_tmp
   real, dimension(2):: q_fit, T_fit
!
   td=nd
   T_fit(1)=T_surf-1.0
   T_fit(2)=T_surf+1.0
   do i=1,2
      e0=2.1718E8*EXP(-4157.0/(T_fit(i)+239.09))
      dbt_tmp = dbt(ncell)
      if (dbt_tmp .lt. 0.0) dbt_tmp = 0.0
!      rb=pf*(dbt(ncell)-T_fit(i))
      rb=pf*(dbt_tmp-T_fit(i))
      lvp=597.0-0.57*T_fit(i)
      q_evap0=1000.*lvp*evap_coeff*wind(ncell)
      q_conv=rb*q_evap0
      q_evap=q_evap0*(e0-ea(ncell))
      !
      if(q_evap.lt.0.0) then 
         q_evap=0.0
         q_conv = -0.65*rb*q_evap0
      end if
      q_ws=6.693E-2+1.471E-3*T_fit(i)

      q_fit(i)=q_ns(ncell)+q_na(ncell)-q_ws-q_evap+q_conv
   end do
!
!     q=AT+B
>>>>>>> origin
!
!     Linear fit over the range of 2.0 deg C.
!     These results can be used to estimate the "equilibrium" 
!     temperature and linear rate constant.
!
   A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
   q_surf=0.5*(q_fit(1)+q_fit(2))
   B=(q_surf/A)-(T_fit(1)+T_fit(2))/2.
!
!     ******************************************************
!               Return to Subroutine RIVMOD
!     ******************************************************
!
END Subroutine Energy
