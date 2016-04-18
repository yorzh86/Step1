module output_mod
	use kinds_mod
	use units_mod
	use system_mod
	use integrate_mod
	implicit none
	
contains

	subroutine showResults
		integer::i
		print *,
		write(*,*) 'Writing time[ps], timestep and temperatures for regions:'
		do i=0, N_steps
			write(*,'(1X, 1F8.2, 1F10.0, 10F15.8 )') regions(i)%temps
		end do

		print *,
		write(*,*) 'Writing time[ps], timestep and kinetic energies before swap:'
		
		!how we store energies:
		do i=0, N_steps
			if(.not. allocated(regions(i)%energies)) cycle
			write(*, '(1X, 1F8.2, 1F10.0, 2F17.10)') regions(i)%energies
		end do
		
	end subroutine showResults

	subroutine writeStepXYZ(iou)
		integer,intent(in)::iou
		
		character(32)::buf
		integer::i,k
		
		write(buf,*) size(atoms)
		write(iou,'(1A)') trim(adjustl(buf))
		write(iou,'(1A,1I10)') 'Atoms. Timestep: ',ts
		do k=1,size(atoms)
			write(iou,'(1I4, 3F5.1,6F13.9)') atoms(k)%atom_id, &
				& [(convert(atoms(k)%r(i),'m','A'),i=1,3)], &
				& [(convert(atoms(k)%v(i),'m/s', 'A/ps'), i=1,3)], &
				& [(convert(atoms(k)%f(i), 'N', 'eV/A'), i=1,3)]            
		end do
	end subroutine writeStepXYZ
	
	subroutine writeStepThermo (k, iou_temps)
		integer, intent(in)::k, iou_temps 
		real(wp)::t
			   
		if (k==0) then
			write(iou_temps, '(1X, 1A8, 2X, 1A8, 2X, 1A20, 1I3, 1A9 )')'Time[ps],', & 
				& 'TimeStep[ms],', 'Temperatures[K] for:', N_slabs, ' regions.'
			write(iou_temps,*)
		end if
		
		t = k*1E-2_wp

		write(iou_temps,'(1X, 1F8.2, 1F10.0, 10F15.8)') regions(k)%temps

	end subroutine writeStepThermo
	
	subroutine writeStepEnergies(k, iou_energies)
		integer, intent(in)::k, iou_energies
		real(wp)::t
		
		t = convert((k*dt),'s','ps')
		
		if (k==0) then      
			write(iou_energies,'(1X, 1A8, 2X, 1A8, 2X, 1A37)') 'Time[ps],', & 
				& 'TimeStep[ms],', 'Kinetic energy[eV] of swapped atoms.'
			write(iou_energies,*)
		end if
		
		write(iou_energies,'(1X, 1F7.2, 1I7, 2F15.8)') t,  k, & 
			& convert(KEi(hot), 'J','eV'), &
			& convert(KEi(cold), 'J', 'eV')
			
	end subroutine writeStepEnergies
	
	subroutine writeBasicInfo ()
		integer::i

		write(*,'(/, 1X,1A12, T29, 1A1, 2(1F4.1,", "), 1F5.1, 1A1)')'Box size[A]:', &
			&  '[',[(convert(box(i), 'm','A'),i=1,3)],']'
		write(*,'(1X, 1A15, T30, 1F5.1)') 'Temperature[K]:', T0
		write(*,'(1X, 1A16, T30, 1I5)') 'Number of atoms:', size(atoms)
		write(*,'(1X, 1A16, T30, 1I6)') 'Number of steps:', N_steps
		write(*,'(1X, 1A17, T30, 1I3,/)') 'Skip_swap factor:', skip_swap

		
	end subroutine writeBasicInfo



end module output_mod
