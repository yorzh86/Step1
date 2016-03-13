module output_mod
    use kinds_mod
    use units_mod
    use system_mod
    use integrate_mod
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
            write(iou,'(1I4, 3F5.1,6F13.9)') atoms(k)%atom_id, &
                & [(convert(atoms(k)%r(i),'m','A'),i=1,3)], &
                & [(convert(atoms(k)%v(i),'m/s', 'A/ps'), i=1,3)], &
                & [(convert(atoms(k)%f(i), 'N', 'eV/A'), i=1,3)]            
        end do
    end subroutine writeStepXYZ
    
    
    subroutine writeStepThermo (k, iou_temps, iou_energies)
        integer, intent(in)::k, iou_temps, iou_energies
        
        integer,dimension(:), allocatable::l
        integer:: h, c, j
        real(wp):: abc, t
        
        abc = real(latM(3)*lattice_const/N_slabs, wp)

        do j=1, N_slabs
            l = regionList(j*abc - abc, j*abc)
            regions(k+1)%temps(j) = listTemp(l)
            if (j==1) h = selectHot(l)
            if (j==5) c = selectCold(l)
        end do
        t = k*1E-2_wp
        write(iou_temps,'(1X, 1F7.2, 1I7, 10F15.8)') t, k, regions(k+1)%temps
        
        if (mod(k,10)==0) then
            write(iou_energies,'(1X, 1F7.2, 1I7, 2F15.8)') t,  k, & 
                & convert(Ei(h), 'J','eV'), &
                & convert(Ei(c), 'J', 'eV')
            call swapAtoms(h,c)
        end if
        
    end subroutine writeStepThermo

end module output_mod
