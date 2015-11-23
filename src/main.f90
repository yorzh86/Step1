! Draw a bounding box in pymol
!
! http://www.pymolwiki.org/index.php/DrawBoundingBox

program main_prg
	use kinds_mod
	use system_mod
	use integrate_mod
	use output_mod
	implicit none
	
	integer,parameter::Ns = 30000
	integer,parameter::skip = 10
	
	real(wp)::T0 = 40.0_wp
	!! initial temperature [K]
	real(wp)::dt = 10.0E-15_wp
	!! timestep [seconds]
	real(wp)::latticeConstant = 5.26E-10_wp
	integer::iou
	
	call setupSim()
	call runSim()
	
contains

	subroutine setupSim
		open(file='out.xyz',newunit=iou)
		enableLennardJones = .true.
		call setThermostat(.true.,T0,100.0_wp*dt)
		call buildSystem(latticeConstant,[5,5,5],T0)
		
		call doBox()
		call writeStepXYZ(iou)
	end subroutine setupSim

	subroutine runSim
		integer::i,k
		write(*,*) "   k   Temperature[K]    KE[units]     PE[units]	   Heat flux[units]"
		do k=0,Ns
			!call velocityVerlet(dt)
			call leapFrog(dt)
			call doBox()
			if(k==Ns/2) call setThermostat(.false.)
						
			if(mod(k,skip)==0) then
				call writeStepXYZ(iou)
				write(*,'(1I5,10EN15.3)') k,temperature(),KE(),PE(), heatflux()
			end if
			if(mod(k,1000)==0) then
				write(*,*) 'Average Neighbors: ', averageNeighbors()
			end if
			
			if(mod(k,20)==0) then
				do i=1,size(atoms)
					call updateNeighbors(i)
				end do
			end if
			
		end do
		
		close(iou)
	end subroutine runSim

end program main_prg 

