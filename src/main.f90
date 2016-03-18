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
	integer::iou_lammps=1
		!! I/O unit to read lammps dump file
	
	call setupSim()
	!call runSim()
	call endSim()
	
contains

	subroutine setupSim
		integer::i, j, k, atom_id
		real(wp)::simple
		real(wp), dimension(3)::posit, velocity, force
		open(file='mark1.xyz',newunit=iou_xyz)
		open(file='mark1.thermo',newunit=iou_thermo)
		!open(file='../lammps/lammps.all', status= 'old', unit= iou_lammps)
				
		simple = 2.0
		enableLennardJones = .true.
		
		call initialize_parameters()
		!! Initialize E0, S0, cutoff, N-steps, etc (settings.f90)
		call setThermostat(.false.,T0,10.0_wp*dt)
		call setBarostat(.false.,P0, 5.0E10_wp*dt)
		do k=2,3
			if (allocated(types)) deallocate(types)
			if (allocated(atoms)) deallocate(atoms)
			call SimpleSystem (convert(simple,'A','m'), k, T0)
			call printSimpleSystem()
		end do
		
		call doBox()
		call writeLammpsData('Ar.data')
		call writeLammpsVars('Ar.vars')
		
		!call buildSystem(convert(lattice_const,'A','m'),[2,2,2],T0)
		!rewind(iou_lammps)

	end subroutine setupSim
		
	subroutine runSim
		integer::k
		do k=0,N_steps
			!call integrateLammps(dt)
			if(mod(k,skip_thermo)==0) call thermoReport(k)
			if(mod(k,skip_dump  )==0) call writeStepXYZ(iou_xyz)
			if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
			call velocityVerlet(dt)
			
			call doBox()
 			if(k==N_steps/2) call setThermostat(.false.)

		end do
	end subroutine runSim

	subroutine endSim
		close(iou_xyz)
		close(iou_thermo)
		!close(iou_lammps)
	end subroutine endSim
	
	subroutine printSimpleSystem()
		integer::i,j
		write(*,*)
		write(*,*) 'Number of atoms:', size(atoms)
		write(*,*)'Box size[A]:',box*1E10_wp
		write(*,*)
		
		write(*,'(1X, 1A5, 6A25, 1A13 )')  'id', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz', 'fnorm'
		
		do i=1, size(atoms)
			write(stdout, '(1X, 1I5, 6E25.15, 1F13.6)') atoms(i)%atom_id, &
				& [(convert(atoms(i)%r(j),'m','A'),j=1,3)], &
				& [(convert(atoms(i)%f(j), 'N', 'eV/A'), j=1,3)], &
				& convert(fnorm(),'N','eV/A')
		end do
	end subroutine printSimpleSystem

	subroutine thermoReport(k)
		integer,intent(in)::k
		integer,save::c = 0
		!character(128)::buf
		real(wp), dimension(3)::o
		integer::i,j
		
		o = heatflux()
		
		if (mod(c,50)==0) then
			write(*,*)
			write(stdout,'(1X,1A17,1A5, 3F6.2)')"\x1B[96mSystem size:","\x1B[97m", [(convert(box(i), 'm', 'A'), i=1,3)]
			write(stdout,'(1X,1A21, 1A5, 1I5)') "\x1B[96mNumber of atoms:","\x1B[97m", size(atoms)
			write(stdout,'(1X,1A29, 1A5, 1I3)') "\x1B[96mAv. number of neighbors:", "\x1B[97m", nint(averageNeighbors())
			write(*,*)			
			write(*, '(1X, 1A5, 1A9, 6A13, 1A22)') 'Step', 'Temp', 'KE()', 'PE()', 'TotEng', 'Jx', 'Jy', 'Jz', 'Fnorm'
		end if
		
		write(stdout,'(1X,1I3,1F11.3,3F13.6,3ES17.6, 1F13.6)') k, temperature(), &
			& convert(KE(),'J','eV'), &
			& convert(PE(),'J','eV'), &
			& convert(E(),'J','eV'),  &
			& (convert(o(i),'W/m2','eV/ps/A2')/product(box),i=1,3), &
			& convert(fnorm(),'N','eV/A')
		c = c+1
	end subroutine thermoReport

end program main_prg 
