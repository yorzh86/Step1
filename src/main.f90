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
	use properties_mod
	use autodiff_mod
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
		
		write(*,'(/,1X, 1A22,/ )') '\x1B[32;1mFINISHED!!!\x1B[0m'
	end subroutine runSim
	
	subroutine endSim
		!call showResults()
		call thermalConductivity()

		call doMessage(5, "Check grad calculation properties.f90 line75", [stdout])
		call doMessage(3, "Re-do updateAllLists in system, check commit Meeting", [stdout])
		call doMessage(3, "Check lammps script with system relaxation, and correct py-script", [stdout])
		call doMessage(5, "Add to code calc dkSI/dEpsilon", [stdout])
		close(iou_xyz)
		close(iou_temps)
		close(iou_energies)
		close(iou_log)
	end subroutine endSim

end program main_prg 
