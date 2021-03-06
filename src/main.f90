!! Draw a bounding box in pymol
!http://www.pymolwiki.org/index.php/DrawBoundingBox

!! Using Advanced metal units

!mass = grams/mole 
!distance = Angstroms 
!time = picoseconds 
!---
!energy = u*A^2/ps^2 (instead of eV)
!---
!velocity = Angstroms/picosecond +
!---
!force = u*A^2/ps^2/Angstrom= u*A/ps^2 (instead of eV/Angstrom)
!---
!torque = u*A^2/ps^2 (instead of eV)
!---
!temperature = Kelvin 
!pressure = bars 
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
	!! I/O unit to write energies before swap
	integer::iou_penergies
	integer::iou_totenergies
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
		open(file='step1.Totenergies',newunit=iou_totenergies)
		open(file='step1.log',newunit=iou_log)

		call initialize_parameters()

		enableLennardJones = .true.
		call setThermostat(.true.,T0,10.0_wp*dt) !turn it ON
		call setBarostat(.true.,P0, 5.0E10_wp*dt)!turn it ON
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
			!if(k==N_steps/3) then
			!	call setThermostat(.false.)
			!	call setBarostat(.false.)
			!end if
			
			call rnem(k)
			
			if(mod(k,skip_swap)==0)     call swapAtoms(k)
			if(mod(k,skip_thermo)==0)   call writeStepThermo(k, iou_temps)
			if(mod(k,skip_swap)==0)     call writeStepEnergies(k,iou_energies, iou_penergies, iou_totenergies)
 			!if(mod(k,skip_dump)==0)     call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllLists()
			
			call velocityVerlet(dt)
			call doBox()
			
		end do
		
		write(*,'(/,1X, 1A22,/ )') '\x1B[32;1mFINISHED!!!\x1B[0m'
	end subroutine runSim
	
	subroutine endSim
		
 		!call showResults()  !output.f90
		!call thermalConductivity() !calculates k
		!call specificHeat() !calculates Cp

		!print *, "below should be same number:"
		!write (*,*) convert(1.0_wp, 'eV', 'uA2/ps2')
		!write (*,*) convert(1.6022E-19_wp, 'J', 'uA2/ps2')
		!write (*,*) convert (1.0_wp, 'eV', 'J')
		
		!write(*,*) "u*A2/ps2 /K"
		!write(*,*) "----------------------"
		!write(*,*) "J/K = kg * m2/s2 /K"


		print *, 
		call doMessage(0, "test message, Check output.f90 doMessage subroutine.", [stdout])
		
		close(iou_xyz)
		close(iou_temps)
		close(iou_energies)
		close(iou_penergies)
		close(iou_totenergies)
		close(iou_log)
	end subroutine endSim

end program main_prg
