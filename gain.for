c     *****************************************************************
c     Program 'gain.for' **********************************************
c     Calculates propagation of sound from a spherical transducer. 
c     *****************************************************************
	implicit real*8 (a-h,o-z)
	real*8 u(0:50000),v(0:50000),rho(0:49999),pressure(0:49999)
	real*8 r(0:50000),volume(0:50000),volume_initial(0:50000)
      pi=2*dacos(0.0d0)
c     ******************************************************************
c     Open file for input. *********************************************
c     ******************************************************************
	open(1,file='gain.inp',status='old')
      read(1,*) frequency
	omega=2*pi*frequency
      read(1,*) rho_initial ! initial density of liquid.  
	read(1,*) r_step
	read(1,*) t_step  ! Time step. 
	read(1,*) t_save  ! Time between saving results as a funciton of r. 
	rho_0=0.14513
	i_t_save=nint(t_save/t_step)
	read(1,*) t_stop  ! Time to stop. 
	i_t_stop=nint(t_stop/t_step)
	read(1,*) v_wall_0	!  Amplitude of velocity at the wall. 
	u_wall_0=v_wall_0/omega  ! Amplitude of displacement at wall. 
	r_max=0.8	   ! Inner radius of transducer. 
	n_r_max=nint(r_max/r_step)
	r_step=r_max/n_r_max
c     ******************************************************************
c     Open file for pressure as a function of r. ***********************
c     ******************************************************************
	open(2,file='gain.dat',status='unknown')
c     ******************************************************************
c     Open file for pressure at center as a function of time. **********
c     ******************************************************************
	open(3,file='gain_focus.dat',status='unknown')
c     ******************************************************************
c     Test subroutine for pressure as a function of density. ***********
c     ******************************************************************
      open(4,file='pressure.dat',status='unknown')
      i_max=600
      do i=0,600
	 rho_test=0.145+0.0001*i
	 call subroutine_rho_pressure(rho_test,pressure_test)
	 write(4,130) rho_test,pressure_test
130    format(f12.5,f12.1)
      end do
	close(4)
c     ******************************************************************
c     Find pressure for bulk density. **********************************
c     ******************************************************************
	call subroutine_rho_pressure(rho_initial,pressure_initial)
      write(*,131) rho_initial,pressure_initial/1.0d6
131   format(f12.5,f12.3)
c     ******************************************************************
c     Set values for r. ************************************************
c     ******************************************************************
	do n_r=0,n_r_max
	 r(n_r)=n_r*r_step
	end do
c     ******************************************************************
c     Set initial values of u and v. ***********************************
c     ******************************************************************
	do n_r=0,n_r_max
	 u(n_r)=0.0d0
	 v(n_r)=0.0d0
	end do
c     ******************************************************************
c     Set initial values of volume_initial(n_r) between r=n_r*d_r 
c     and r=(n_r+1)*d_r. ***********************************************
c     ******************************************************************
	do n_r=0,n_r_max-1 
       volume_initial(n_r)=(4*pi/3)*(r(n_r+1)**3-r(n_r)**3)
	 rho(n_r)=rho_initial
	end do
c     ******************************************************************
c     Begin loop over time. ********************************************
c     ******************************************************************
	pressure_min=0.0d0
      do i_t=0,i_t_stop
	 t=i_t*t_step
c      *****************************************************************
c      Update density. *************************************************
c      *****************************************************************
	 do n_r=0,n_r_max-1 
        volume(n_r)=(4*pi/3)*
     &  ((r(n_r+1)+u(n_r+1))**3-(r(n_r)+u(n_r))**3)
        rho(n_r)=rho_initial*volume_initial(n_r)/volume(n_r)
        call subroutine_rho_pressure(rho(n_r),pressure(n_r))
	 end do
c      *****************************************************************
c      Update velocity. ************************************************
c      *****************************************************************
	 v(0)=0.0d0
       do n_r=1,n_r_max-1 
	  v(n_r)=v(n_r)-(pressure(n_r)-pressure(n_r-1))*2*t_step/
     &  ((rho(n_r)+rho(n_r-1))*r_step)
	 end do
	 call subroutine_wall(t,omega,u_wall_0,
     & u_wall,v_wall)  
  	 v(n_r_max)=v_wall
c      *****************************************************************
c      Update displacement. ********************************************
c      *****************************************************************
	 u(0)=0.0d0
       do n_r=1,n_r_max-1 
	  u(n_r)=u(n_r)+v(n_r)*t_step
	 end do
	 call subroutine_wall(t,omega,u_wall_0,u_wall,v_wall)  
  	 u(n_r_max)=u_wall
c      *****************************************************************
c      Smooth displacement. ********************************************
c      *****************************************************************
       f=0.0001*t_step/r_step
c       do n_r=1,n_r_max-1
c	  u_smooth(n_r)=f*u(n_r-1)+(1-2*f)*u(n_r)+f*u(n_r+1)
c       end do     
c       do n_r=1,n_r_max-1
c	  u(n_r)=u_smooth(n_r)
c	 end do  
c    	 *****************************************************************
c      Save  results as a function of r. *******************************
c      *****************************************************************
       if (mod(i_t,i_t_save) .eq. 0) then 
	  do n_r=0,n_r_max-1
	   write(2,100) (n_r+0.5)*r_step,(u(n_r)+u(n_r+1))/2,
     &   (v(n_r)+v(n_r+1))/2,rho(n_r),pressure(n_r)/1.0d6
100  	   format(f12.5,f10.6,f10.6,f10.6,f10.6)
	  end do
	  write(2,100)
	 endif  
c    	 *****************************************************************
c      Save  results as a function of t. *******************************
c      *****************************************************************
       if (mod(i_t,i_t_stop/1000) .eq. 0) then 
	  write(3,110) t/1.0d-6,pressure(0)/1.0d6
	  !write(*,110) t/1.0d-6,pressure(0)/1.0d6
110  	  format(f12.5,f12.6)
	 endif  
c      *****************************************************************
c      Find minimum pressure. ******************************************
c      *****************************************************************
       pressure_min=dmin1(pressure(0),pressure_min) 
c      *****************************************************************
c      Stop. ***********************************************************
c      *****************************************************************
       if (i_t .gt. i_t_stop) goto 200
      end do
200   continue
      write(*,201) pressure_min/1.0d6
201   format('  Pressure_min=',f7.3) 
	stop
	end
c     ******************************************************************
c     ******************************************************************
      subroutine subroutine_rho_pressure(rho,pressure)
c     Calculates the pressure from the density. ************************
c     ******************************************************************
c     ******************************************************************
	implicit real*8 (a-h,o-z)
	c_0=2.383d4
	rho_0=0.14513
	u_0=2.84
	w_0=0.194
	der_1=c_0**2
	der_2=2*c_0**2*u_0/rho_0
	der_3=2*c_0**2*u_0**2/rho_0**2+2*c_0**2*w_0/rho_0**2
     	pressure=(rho-rho_0)*der_1+(rho-rho_0)**2*der_2/2+
     &(rho-rho_0)**3*der_3/6	      
	return 
	end
c     ******************************************************************
c     ******************************************************************
      subroutine subroutine_wall(t,omega,u_wall_0,u_wall,v_wall)
c     Calculates the displacement and velocity of the transducer wall at 
c     time t. The velocity must be consistent with the displacement.  
c     ******************************************************************
c     ******************************************************************
	implicit real*8 (a-h,o-z)
	tau=10.0d-6   ! Time for transducer amplitude to build up.  
      u_wall=u_wall_0*dsin(omega*t)*(1-dexp(-t/tau))
      v_wall=u_wall_0*omega*dcos(omega*t)*(1-dexp(-t/tau))+
     &dsin(omega*t)*dexp(-t/tau)/tau
	return 
	end


