!! Draw a bounding box in pymol
!http://www.pymolwiki.org/index.php/DrawBoundingBox

!! Back to metal units
!mass = grams/mole +
!distance = Angstroms +
!time = picoseconds +
!energy = eV +
!velocity = Angstroms/picosecond +
!force = eV/Angstrom +
!torque = eV
!temperature = Kelvin
!pressure = bars +
!dynamic viscosity = Poise
!charge = multiple of electron charge (1.0 is a proton)
!dipole = charge*Angstroms
!electric field = volts/Angstrom
!density = gram/cm^dim

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
	integer::iou_penergies
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
		open(file='step1.Penergies',newunit=iou_penergies)
		open(file='step1.log',newunit=iou_log)

		call initialize_parameters()

		enableLennardJones = .true.
		call setThermostat(.false.,T0,10.0_wp*dt) !turn it ON
		call setBarostat(.false.,P0, 5.0E10_wp*dt)!turn it ON
		call buildSystem(lattice_const,latM,T0)
		call doBox()
		call writeBasicInfo()
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
		print *, "system setup successful"
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
			!if(mod(k,skip_thermo)==0)   call writeStepThermo(k, iou_temps)
			!if(mod(k,skip_swap)==0)     call writeStepEnergies(k,iou_energies, iou_penergies)
 			!if(mod(k,skip_dump)==0)     call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllLists()
			
			
			!if(mod(k,1)==0) call test_diff(k) !output.f90 prints 
			call velocityVerlet(dt)
			call doBox()
		end do
		
		write(*,'(/,1X, 1A22,/ )') '\x1B[32;1mFINISHED!!!\x1B[0m'
	end subroutine runSim
	
	subroutine endSim
 		!call showResults()  !output.f90
		!call thermalConductivity() !calculates k
		
		call doMessage(0, "test message, Check output.f90 doMessage subroutine.", [stdout])
		print *, 
		call doMessage(0, "check log files: step1.xyz, step1.energies, step1.penergies.", [stdout])
		call doMessage(0, "all log files are in converted units for build: 'build'.", [stdout])
		

		close(iou_xyz)
		close(iou_temps)
		close(iou_energies)
		close(iou_penergies)
		close(iou_log)
	end subroutine endSim

end program main_prg
