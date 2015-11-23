module output_mod
	use kinds_mod
	use system_mod
	implicit none
	
contains

	subroutine writeStepXYZ(iou)
		integer,intent(in)::iou
		
		character(32)::buf
		integer::k
		
		write(buf,*) size(atoms)
		write(iou,'(1A)') trim(adjustl(buf))
		write(iou,'(1A,1I10)') 'Atoms. Timestep: ',ts
		do k=1,size(atoms)
			write(iou,'(1A2,1X,3F25.15)') types(atoms(k)%t)%atom_name,atoms(k)%r*1E10_wp
			!! multiply by 1E10 to secure the correct work of pymol
		end do
	end subroutine writeStepXYZ

end module output_mod 
