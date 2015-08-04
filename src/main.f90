program main_prg
	use kinds_mod
	use system_mod
	use integrate_mod
	use output_mod
	implicit none
	
	integer,parameter::Ns = 2000
	integer,parameter::skip = 1
	
	real(wp),dimension(Ns/skip,6)::plotData
	real(wp)::T0 = 40.0_wp
	real(wp)::dt = 5.0_wp
	integer::iou
	
	call setupSim()
	call runSim()
	call doPlots()
	
contains

	subroutine setupSim
		open(file='out.xyz',newunit=iou)
		
		enableLennardJones = .true.
		call setThermostat(.true.,T0,100.0_wp)
		
		call buildSystem(5.260_wp,15,T0)
		!(lattice parameter, box edge, temperature)
		call doBox()
		call writeStepXYZ(iou)
	end subroutine setupSim

	subroutine runSim
		integer::k
		
		do k=1,Ns
			call velocityVerlet(dt)
			call doBox()
			if(k==1000) call setThermostat(.false.)
			
			if(mod(k,skip)==0) then
				call writeStepXYZ(iou)
				write(*,'(1I5,10ES15.3)') k,temperature(),PE(),KE(),heatflux()
				plotData(k/skip,:) = [t,temperature(),PE(),KE(),heatflux()]
			end if
		end do
		
		close(iou)
	end subroutine runSim

	subroutine doPlots
		use plplotlib_mod
		
		real(wp),dimension(:),allocatable::t,Tp,KE,PE,Jx,Jy
		real(wp),dimension(:),allocatable::JJx,JJy,IJJ
		
		t  = plotData(:,1)
		Tp = plotData(:,2)
		PE = plotData(:,3)
		KE = plotData(:,4)
		Jx = plotData(:,5)
		Jy = plotData(:,6)
		
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
		call xylim(mixval(t),mixval(IJJ))
		call plot(t,IJJ,lineColor='b',lineWidth=2.0_wp)
		call ticks()
		call labels('Integration Variable #fi#gt#fn [ps]','Cumulative ACF Inegral [TODO]','')
		
		call show()
	end subroutine doPlots

	function autocorrelate(A) result(o)
		real(wp),dimension(:),intent(in)::A
		real(wp),dimension(:),allocatable::o
		
		integer::i,j
		
		allocate(o(size(A)))
		
		do i=1,size(A)
			o(i) = 0.0
			do j=i,size(A)
				o(i) = o(i)+A(j)*A(j-i)
			end do
			o(i) = o(i)/real(size(A)-1,wp)
		end do
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

