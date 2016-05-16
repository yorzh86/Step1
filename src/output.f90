module output_mod
	use kinds_mod
	use units_mod
	use system_mod
	use integrate_mod
	implicit none
	
contains

	subroutine showResults
		integer::i, j

		print *,
		write(*,*) 'Writing temperatures for regions:'
		
		do i=0, N_steps
			write(*,'(1X, 10F15.8)') (regions(j)%temps(i)%x, j=1,10)
		end do
		
		write(*,*)
		write(*,*) 'Writing KEi of two atoms:'
		do i=0, N_steps
			!if(.not. allocated( (regions(j)%energies(i),j=1,2 ))) cycle !!!!!!!!!!does not work
			write(*, *) (regions(j)%energies(i)%x, j=1,2)
		end do
		write(*,*)

	end subroutine showResults

	subroutine writeStepXYZ(iou)
		integer,intent(in)::iou
		
		character(32)::buf
		integer::i,k
		
		write(buf,*) size(atoms)
		write(iou,'(1A)') trim(adjustl(buf))
		write(iou,'(1A,1I10)') 'Atoms. Timestep: ',ts
		do k=1,size(atoms)
			write(iou,'(1I4, 3F5.1,6E25.15)') atoms(k)%atom_id, &
				& [(convert(atoms(k)%r(i),'m','A'),i=1,3)], &
				& [(convert(atoms(k)%v(i),'m/s', 'A/ps'), i=1,3)], &
				& [(convert(atoms(k)%f(i), 'N', 'eV/A'), i=1,3)]            
		end do
	end subroutine writeStepXYZ
	
	subroutine writeStepThermo (k, iou_temps)
		integer, intent(in)::k, iou_temps 
		real(wp)::t
		integer::j,i
			   
		if (k==0) then
			write(iou_temps, '(1X, 1A8, 2X, 1A8, 2X, 1A20, 1I3, 1A9 )')'Time[ps],', & 
				& 'TimeStep[ms],', 'Temperatures[K] for:', N_slabs, ' regions.'
			write(iou_temps,*)
		end if
		
		t = convert(real(k,wp)*dt,'s','ps')
		
		write(iou_temps,'(1X, 1F8.2, 1F10.0, 10F15.8)') t, real(ts, wp), (regions(j)%temps(k)%x, j=1,10)

	end subroutine writeStepThermo
	
	subroutine writeStepEnergies(k, iou_energies)
		integer, intent(in)::k, iou_energies
		real(wp)::t
		
		t = convert(real(k,wp)*dt,'s','ps')
		
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

		write(*,'(1X,1A26,T35,1A1,2(1F4.1,", "),1F5.1,1A1)')'\x1B[37;1mBox size[A]:\x1B[33;1m:', &
			&  '[',[(convert(box(i), 'm','A'),i=1,3)],']'
		write(*,'(1X, 1A29, T35, 1A6)') '\x1B[37;1mTemperature[K]:\x1B[33;1m',&
			& adjustl(real2char(real(T0)))
		write(*,'(1X, 1A30, T35, 1A6)')  '\x1B[37;1mNumber of atoms:\x1B[33;1m',& 
			& adjustl(int2char(size(atoms)))
		write(*,'(1X, 1A30, T35, 1A6)')  '\x1B[37;1mNumber of steps:\x1B[33;1m',&
			& adjustl(trim(int2char(N_steps)))
		write(*,'(1X, 1A31, T35, 1A6)')'\x1B[37;1mUpdate neighbors:\x1B[33;1m',&
			& adjustl(int2char(skip_neighbor))
		write(*,'(1X, 1A31, T35, 1A6,/)')'\x1B[37;1mSkip_swap factor:\x1B[33;1m',&
			& adjustl(int2char(skip_swap))
		print *, '\x1B[0m'

	end subroutine writeBasicInfo


	subroutine doMessage(priority, msg, ious)
		integer, intent(in)::priority
		character(len=*), intent(in)::msg
		integer, dimension(:),intent(in)::ious
		character(12) :: time
		integer,dimension(8) :: values
		integer::k

		character(:), allocatable::temp1
		character(:), allocatable::color1
		character(:), allocatable::color2
		
		call date_and_time(VALUES=values)
		time = '<'//int2char(values(5))//':'//int2char(values(6))//':'//int2char(values(7))//'>'

		select case (priority)
			case(1)
				temp1='Error: '
				color1 = '\x1B[31;1m'
				color2 = '\x1B[37;1m'
			case(2)
				temp1='Warning: '
				color1 = '\x1B[33;1m'
				color2 = '\x1B[37;1m'
			case(3)
				temp1='Debug: '
				color1 = '\x1B[35;1m'
				color2 = '\x1B[37;1m'
			case default
				temp1='Attention: '
				color1 = '\x1B[36;1m'
				color2 = '\x1B[37;1m'
		end select
		
		do k=1,size(ious)
			if(ious(k)==stdout) then
				write(*,*) color1//time, color1//temp1//color2, msg
				print *, '\x1B[0m'
			else
				write(ious(k),*) time, temp1, msg
			end if
		end do
	end subroutine doMessage

end module output_mod
