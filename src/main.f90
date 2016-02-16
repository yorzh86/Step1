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
	implicit none
	
	integer::iou_xyz
		!! I/O unit for xyz file output
	integer::iou_thermo
		!! I/O unit for thermo report output
	
	call setupSim()
	call runSim()
	call endSim()
	
contains

	subroutine setupSim
		integer::i, j, k, atom_id
		real(wp), dimension(3)::posit, velocity, force
		open(file='mark1.xyz',newunit=iou_xyz)
		open(file='mark1.thermo',newunit=iou_thermo)
		open(1, file='../lammps/lammps.all', status= 'old')
		
		
		
		enableLennardJones = .true.
		
		call initialize_parameters()
		!! Initialize E0, S0, cutoff, N-steps, etc (settings.f90)
		call setThermostat(.false.,T0,10.0_wp*dt)
		call setBarostat(.false.,P0, 5.0E10_wp*dt)
		call buildSystem(convert(lattice_const,'A','m'),[2,2,2],T0)
		rewind(1)
		call doBox()
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
		
	end subroutine setupSim
		
	subroutine runSim
		integer::k
		open(1, file='../lammps/lammps.all', status= 'old')
		do k=0,N_steps
			call integrateLammps(dt)
			if(mod(k,skip_thermo)==0) call thermoReport(k)
			if(mod(k,skip_dump  )==0) call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
			
			!call velocityVerlet(dt)
			call doBox()
 			if(k==N_steps/2) call setThermostat(.false.)

		end do
	end subroutine runSim

	subroutine endSim
		close(iou_xyz)
		close(iou_thermo)
		close(1)
	end subroutine endSim

	subroutine thermoReport(k)
		integer,intent(in)::k
		integer,save::c = 0
		!character(128)::buf
		real(wp), dimension(3)::o
		integer::i	
		
		o = heatflux()
		
		if (mod(c,50)==0) then
!			write(*,*)
!			write(stdout,'(1X,1A17,1A5, 3F6.2)')"\x1B[96mSystem size:","\x1B[97m", [(convert(box(i), 'm', 'A'), i=1,3)]
!			write(stdout,'(1X,1A21, 1A5, 1I5)') "\x1B[96mNumber of atoms:","\x1B[97m", size(atoms)
!			write(stdout,'(1X,1A29, 1A5, 1I3)') "\x1B[96mAv. number of neighbors:", "\x1B[97m", nint(averageNeighbors())
!			write(*,*)			
!			write (stdout, '(1X, 1A5, 1A4, 1A8, 1A10, 1A11, 2A18, 1A4)') "\x1B[93m", 'k[#]','temperature', &
!				& 'tEnergy', 'Jx','Jy','Jz',  "\x1B[0m"
		write(*, '(1X, 1A5, 1A10, 6A13)') 'Step', 'Temp', 'KE()', 'PE()', 'TotEng', 'Jx', 'Jy', 'Jz'
		end if

		write(stdout,'(1X,1I3,1F11.3,3F13.6,3ES17.6, 1F13.6)') k, temperature(), &
			& convert(KE(),'J','eV'), &
			& convert(PE(),'J','eV'), &
			& convert(E(),'J','eV'),  &
			& (convert(o(i),'W/m2','eV/ps/A2')/product(box),i=1,3), &
			& convert(sumdelV(),'N','eV/A')
		
		!write(stdout,'(1X,1I3,1F11.3,3F13.5,3ES17.6)') k, temperature(), E(), KE(), PE(), o(i)/product(box)
		
		c = c+1
	end subroutine thermoReport

end program main_prg 
