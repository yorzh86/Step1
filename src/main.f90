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
	implicit none
	
	integer::iou_xyz
		!! I/O unit for xyz file output
	integer::iou_energies
		!! I/O unit to write energies before swap
	integer::iou_temps
		!! I/O unit to write temperatures
	integer::iou_log
		!! I/O unit for log_file
	real(wp)::p

	call setupSim()
	call runSim()
	call endSim()
	
contains

	subroutine setupSim

		open(file='step1.xyz',newunit=iou_xyz)        
		open(file='step1.temps',newunit=iou_temps)
		open(file='step1.energies',newunit=iou_energies)
		open(file='step1.log',newunit=iou_log)

		call initialize_parameters()

		enableLennardJones = .true.
		call setThermostat(.true.,T0,10.0_wp*dt)
		call setBarostat(.true.,P0, 5.0E10_wp*dt)
		call buildSystem(lattice_const,latM,T0)
		call doBox()
		call writeBasicInfo()
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
	end subroutine setupSim

	subroutine runSim
		integer::k

		do k=0, N_steps
			p = real(k,wp)/real(N_steps, wp)
			call showProgress(' Simulation ongoing', p)
			if(k==N_steps/3) then
				call setThermostat(.false.)
				call setBarostat(.false.)
			end if

			call rnem(k)

			if(mod(k,skip_swap)==0)     call swapAtoms(k)
			if(mod(k,skip_thermo)==0)   call writeStepThermo(k, iou_temps)
			if(mod(k,skip_swap)==0)     call writeStepEnergies(k,iou_energies)
			!if(mod(k,skip_dump)==0)    call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllLists()
			
			call velocityVerlet(dt)
			call doBox()
		end do
		
		write(*,'(/,1X, 1A23,/ )') '\x1B[37;42mFINISHED!!!\x1B[0m'
	end subroutine runSim
	
	subroutine endSim
	integer::i

	do i=1, 4
		call doMessage(i, 'hello world!', .true., iou_log)
	end do
		close(iou_xyz)
		close(iou_temps)
		close(iou_energies)
		close(iou_log)
	end subroutine endSim

end program main_prg 
