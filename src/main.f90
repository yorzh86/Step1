program main_prg
	use kinds_mod
	use system_mod
	use integrate_mod
	use output_mod
	implicit none
	
	integer,parameter::Ns = 10000
	integer,parameter::skip = 10
	
	real(wp),dimension(Ns/skip,6)::plotData
	real(wp)::T0 = 60.0_wp
	!! initial temperature [K]
	real(wp)::dt = 10E-15_wp
	!! timestep [seconds]
	real(wp)::latticeConstant = 5.26E-10_wp
	!real(wp)::latticeConstant = 5.00E-10_wp
	integer::iou
	
	call setupSim()
	call runSim()
! 	call doPlots(plotData)
	
contains

	subroutine setupSim
		open(file='out.xyz',newunit=iou)
		enableLennardJones = .true.
		call setThermostat(.true.,T0,100.0_wp*dt)
		call buildSystem(latticeConstant,8,T0)
		!(lattice parameter, box edge, temperature)
		
		call doBox()
		call writeStepXYZ(iou)
	end subroutine setupSim

	subroutine runSim
		integer::i,k
		write(*,*) box
		write(*,*) "   k   Temperature[K]    KE[units]     PE[units]	   Heat flux[units]"
		do k=0,Ns
			call velocityVerlet(dt)
			!call leapFrog(dt)
			call doBox()
			if(k==Ns/2) call setThermostat(.false.)
			!if(k==5) call setThermostat(.false.)
						
			if(mod(k,skip)==0) then
				call writeStepXYZ(iou)
				write(*,'(1I5,10EN15.3)') k,temperature(),KE(),PE(), heatflux()
! 				plotData(k/skip,:) = [t,temperature(),KE(),PE(),heatflux()]
			end if
			
			if(mod(k,20)==0) then
				do i=1,size(atoms)
					call updateNeighbors(i)
				end do
			end if
			
		end do
		
		close(iou)
	end subroutine runSim

	subroutine doPlots(plotData)
		use plplotlib_mod
		real(wp),dimension(:,:),intent(in)::plotData
		
		real(wp),dimension(:),allocatable::t,Tp,KE,PE,Jx,Jy
		real(wp),dimension(:),allocatable::JJx,JJy,IJJ
		
		integer::N
		
		N = size(plotData,1)
		
		t  = plotData(N/2:,1)
		Tp = plotData(N/2:,2)
		PE = plotData(N/2:,3)
		KE = plotData(N/2:,4)
		Jx = plotData(N/2:,5)
		Jy = plotData(N/2:,6)
		
		N = N/2
		t = t-t(1)
		
		JJx = autocorrelate(Jx)
		JJy = autocorrelate(Jy)
		
		IJJ = cumtrapz( (JJx+JJy)/2.0_wp ,t)
		
		call setup(device='pdfqt',fileName='plot-%n.pdf')
		
		call figure()
		call subplot(1,1,1)
		call xylim(mixval(t),mixval(Tp))
		call plot(t,Tp,lineColor='r',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Temperature #fiT#fn [K]','')
		
		call figure()
		call subplot(2,1,1)
		call xylim(mixval(t),mixval(PE))
		call plot(t,PE,lineColor='b',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Potential Energy #fiPE#fn [eV]','')
		call subplot(2,1,2)
		call xylim(mixval(t),mixval(KE))
		call plot(t,KE,lineColor='r',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Kinetic Energy #fiKE#fn [eV]','')
		
		call figure()
		call subplot(2,1,1)
		call xylim(mixval(t),mixval(Jx))
		call plot(t,Jx,lineColor='b',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Flux #fiJ#dx#u#fn [eV/A#u2#d.ps]','')
		call subplot(2,1,2)
		call xylim(mixval(t),mixval(Jy))
		call plot(t,Jy,lineColor='r',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Flux #fiJ#dy#u#fn [eV/A#u2#d.ps]','')
		
		call figure()
		call subplot(2,1,1)
		call xylim(mixval(t),mixval(JJx))
		call plot(t,JJx,lineColor='b',lineWidth=2.0_wp)
		call ticks()
		call labels('Delay #fi#gd#fn [ps]','Flux ACF #fi<J#dx#u,J#dx#u>#fn [TODO]','')
		call subplot(2,1,2)
		call xylim(mixval(t),mixval(JJy))
		call plot(t,JJy,lineColor='r',lineWidth=2.0_wp)
		call ticks()
		call labels('Delay #fi#gd#fn [ps]','Flux ACF #fi<J#dy#u,J#dy#u>#fn [TODO]','')
		
		
		call figure()
		call subplot(1,1,1)
		call xylim(mixval(t(1:N/2)),mixval(IJJ(1:N/2)))
		call plot(t(1:N/2),IJJ(1:N/2),lineColor='b',lineWidth=2.0_wp)
		call ticks()
		call labels('Integration Variable #fi#gt#fn [ps]','Cumulative ACF Inegral [TODO]','')
		
		call show()
	end subroutine doPlots

	function autocorrelate(A) result(o)
		real(wp),dimension(:),intent(in)::A
		real(wp),dimension(:),allocatable::o
		
		integer::N,i
		
		N = size(A)
		allocate(o(N))
		
		forall(i=1:N) o(i) = sum(A(i:N)*A(1:N-i+1))/real(N-i+1,wp)
	end function autocorrelate

	function cumtrapz(A,t) result(o)
		real(wp),dimension(:),intent(in)::A,t
		real(wp),dimension(:),allocatable::o
		
		integer::k
		allocate(o(size(A)))
		
		o(1) = 0.0_wp
		do k=2,size(A)
			o(k) = o(k-1)+(A(k)+A(k-1))/2.0_wp*(t(k)-t(k-1))
		end do
	end function cumtrapz

end program main_prg 

