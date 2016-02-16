module output_mod
	use kinds_mod
	use units_mod
	use system_mod
	implicit none
	
contains

	subroutine writeStepXYZ(iou)
		integer,intent(in)::iou
		
		character(32)::buf
		integer::i,k
		
		write(buf,*) size(atoms)
		write(iou,'(1A)') trim(adjustl(buf))
		write(iou,'(1A,1I10)') 'Atoms. Timestep: ',ts
		do k=1,size(atoms)
			write(iou,'(1I4,1X,3F13.9)') atoms(k)%atom_id, [(convert(atoms(k)%v(i), 'm/s', 'A/ps'), i=1,3)]
			!write(iou,'(1A2,1X,3F25.15)') types(atoms(k)%t)%atom_name, (atoms(k)%a(i)*types(atoms(k)%t)%m*6.02213665168E26_wp &
			!& * 1.0E10_wp/1.0E12_wp/1.0E12_wp, i=1,3)
			!write(iou,'(1A2,1X,3F25.15, 1F10.5)') types(atoms(k)%t)%atom_name, [(convert(atoms(k)%f(i), 'N', 'eV/A'), i=1,3)]
			
		end do
	end subroutine writeStepXYZ

end module output_mod 
