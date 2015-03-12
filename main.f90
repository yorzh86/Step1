program main_prg
	use kinds_mod
	implicit none
	
	call testOrbit
	
contains

	subroutine testOrbit
		use system_mod
		use integrate_mod
		use output_mod
		
		integer::iou1,iou2
		real(wp)::dt,T0
		
		enableLennardJones = .true.
		
		call buildOrbit(T0)
		dt = T0/100.0_wp
		
		open(file='trajectory.dat',newunit=iou1)
		open(file='trajectory.xyz',newunit=iou2)
		
		call writeOrbit(iou1)
		call writeStepXYZ(iou2)
		do while(t<=T0+dt/2.0_wp)
		
			call velocityVerlet(dt/2.0_wp)
			call leapFrog(dt/2.0_wp)
			call doBox
			
			call writeOrbit(iou1)
			call writeStepXYZ(iou2)
			
		end do
		
		close(iou1)
		close(iou2)
		
		call execute_command_line('./post.py')
	end subroutine testOrbit

end program main_prg 
