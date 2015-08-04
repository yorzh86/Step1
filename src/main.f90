program main_prg
	use kinds_mod
	use system_mod
	use integrate_mod
	use output_mod
	implicit none
	
	integer,parameter::Ns = 1000
	integer,parameter::skip = 5
	
	real(wp),dimension(Ns/skip,6)::plotData
	real(wp)::T0 = 40.0_wp
	real(wp)::dt
	integer::iou
	
	call setupSim()
	call runSim()
	call doPlots()
	
contains

	subroutine setupSim
!~ 		open(file='out.xyz',newunit=iou)
		
		enableLennardJones = .true.
		call setThermostat(.true.,T0,100.0_wp)
		dt = 1.0_wp
		
		call buildSystem(5.260_wp,15,T0)
		!(lattice parameter, box edge, temperature)
		call doBox()
!~ 		call writeStepXYZ(iou)
	end subroutine setupSim

	subroutine runSim
		integer::k
		
		do k=1,Ns
			call velocityVerlet(dt)
			call doBox()
			if(k==10000) call setThermostat(.false.)
			
			if(mod(k,skip)==0) then
!~ 				call writeStepXYZ(iou)
				write(*,'(1I5,10ES15.3)') k,temperature(),PE(),KE(),heatflux()
				plotData(k/skip,:) = [t,temperature(),PE(),KE(),heatflux()]
			end if
		end do
		
		close(iou)
	end subroutine runSim

	subroutine doPlots
		use plplotlib_mod
		
		real(wp),dimension(:),allocatable::t,Tp,KE,PE,Jx,Jy
		
		t  = plotData(:,1)
		Tp = plotData(:,2)
		KE = plotData(:,3)
		PE = plotData(:,4)
		Jx = plotData(:,5)
		Jy = plotData(:,6)
		
		call setup(device='pdfqt',fileName='plot-%n.pdf')
		
		call figure()
		call subplot(1,1,1)
		call xylim(mixval(t),mixval(Tp))
		call plot(t,Tp,lineColor='r',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Temperature #fiT#fn [K]','')
		
		call figure()
		call subplot(1,1,1)
		call xylim(mixval(t),mixval([PE,KE]))
		call plot(t,PE,lineColor='b',lineWidth=2.0_wp)
		call plot(t,KE,lineColor='r',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Energy #fiE#fn [K]','')
		
		call figure()
		call subplot(1,1,1)
		call xylim(mixval(t),mixval([Jx,Jy]))
		call plot(t,Jx,lineColor='b',lineWidth=2.0_wp)
		call plot(t,Jy,lineColor='r',lineWidth=2.0_wp)
		call ticks()
		call labels('Time #fit#fn [ps]','Flux #fiJ#fn [eV/A#u2#d.ps]','')
		
		call show()
	end subroutine doPlots

end program main_prg 

