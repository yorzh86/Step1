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
	

	real(wp)::latticeConstant = 5.26E-10_wp
		!! Lattice spacing between repetitions
	integer::iou_xyz
		!! I/O unit for xyz file output
	integer::iou_thermo
		!! I/O unit for thermo report output
	
	call setupSim()
	call runSim()
	call endSim()
	
contains

	subroutine setupSim
		open(file='out.xyz',newunit=iou_xyz)
		open(file='thermo.log',newunit=iou_thermo)
		
		call initialize()
		
		enableLennardJones = .true.
		call setThermostat(.true.,T0,100.0_wp*dt)
		call buildSystem(latticeConstant,[5,5,5],T0)
		
		call doBox()
		call writeStepXYZ(iou_xyz)
		call writeLammpsData('Ar.data')
		call writeLammpsVars('vars.inc')
	end subroutine setupSim

	subroutine runSim
		integer::k
		
		do k=0,N_steps-1
			
			call velocityVerlet(dt)
			call doBox()
			
! 			if(k==Ns/2) call setThermostat(.false.)
						
			if(mod(k,skip_thermo)==0) call thermoReport(k)
			if(mod(k,skip_dump  )==0) call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
			
		end do
	end subroutine runSim

	subroutine endSim
		close(iou_xyz)
		close(iou_thermo)
	end subroutine endSim

	subroutine thermoReport(k)
		integer,intent(in)::k
		integer,save::c = 0
		character(128)::buf
		
		if(mod(c,20)==0) then
			write(buf,'(1A5,4A15)') 'k [#]','t [ps]','T [K]','KE [eV]','PE [eV]'
			write(stdout,'(2A,1G10.2)') colorize(trim(buf),[5,5,0]),' Nc = ',averageNeighbors()
			write(iou_thermo,'(1A)') '#'//trim(buf)
		end if
		write(stdout,'(1I5,4G15.7)') k,convert(t,'s','ps'),temperature(),convert(KE(),'J','eV'),convert(PE(),'J','eV')
		write(iou_thermo,'(1I5,4G15.3)') k,convert(t,'s','ps'),temperature(),convert(KE(),'J','eV'),convert(PE(),'J','eV')
		
		c = c+1
	end subroutine thermoReport

end program main_prg 

