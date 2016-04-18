! Draw a bounding box in pymol
!
! http://www.pymolwiki.org/index.php/DrawBoundingBox

program main_prg
	use kinds_mod
	use units_mod
	use settings_mod
	use system_mod
	use integrate_mod
	use output_mod
	use lmpIntegrator_mod
    use progress_mod
	implicit none
	
	integer::iou_xyz
		!! I/O unit for xyz file output
	integer::iou_thermo
		!! I/O unit for thermo report output
	integer::iou_energies
		!! I/O unit to write energies before swap
	integer::iou_temps
		!! I/O unit to write temperatures
    real(wp)::p
		
	call setupSim()
	call runSim()
	call endSim()
  call showResults()
	
contains

	subroutine setupSim
		
		open(file='step1.xyz',newunit=iou_xyz)        
		open(file='step1.temps',newunit=iou_temps)
		open(file='step1.energies',newunit=iou_energies)
		
		call initialize_parameters()
					
		enableLennardJones = .true.
		call setThermostat(.true.,T0,10.0_wp*dt)
		call setBarostat(.true.,P0, 5.0E10_wp*dt)
! 		call writeBasicInfo()
		call buildSystem(lattice_const,latM,T0)
		call doBox()
				  
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
	end subroutine setupSim
		
	subroutine runSim
		integer::k, i
		real(wp)::start, finish
		character(128)::buf
		
		
		call cpu_time(start)
		do k=0, N_steps
			p = real(k,wp)/real(N_steps, wp)
			call showProgress1(' Simulation ongoing', p)
			if(k==N_steps/3) then
				call setThermostat(.false.)
				call setBarostat(.false.)
			end if
			
			call rnem(k)
			if(mod(k,skip_swap)==0) 	call swapAtoms(k)
			if(mod(k,skip_thermo)==0)   call writeStepThermo(k, iou_temps)
			if(mod(k,skip_swap)==0)     call writeStepEnergies(k,iou_energies)
			!if(mod(k,skip_dump)==0)    call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
			
			call velocityVerlet(dt)
			call doBox()
		end do
		
		call cpu_time(finish)
		write(*,'(/,1X,1A23, T40, 1A50  )')'FINISHED! Elapsed time:', &
			& writeUsedTime(finish-start)
		!real2time(finish-start)
		
		
	end subroutine runSim
	
	subroutine endSim
		close(iou_xyz)
		close(iou_temps)
		close(iou_energies)
	end subroutine endSim
	
end program main_prg 
