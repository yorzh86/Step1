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
    
    subroutine writeStepThermo (k, iou_temps)
        integer, intent(in)::k, iou_temps 
        integer:: j
        real(wp)::t
               
        if (k==0) then
            write(iou_temps, '(1X, 1A8, 2X, 1A8, 2X, 1A20, 1I3, 1A9 )')'Time[ps],', & 
                & 'TimeStep[ms],', 'Temperatures[K] for:', N_slabs, ' regions.'
            write(iou_temps,*)
        end if
        
        t = k*1E-2_wp

        write(iou_temps,'(1X, 1F7.2, 1I7, 10F15.8)') t, k, regions(k+1)%temps

    end subroutine writeStepThermo
    
    subroutine writeStepEnergies(k, iou_energies, h, c)
        integer, intent(in)::k, iou_energies, h, c
        real(wp)::t
        integer::j
        
        t = k*1E-2_wp
        
        if (k==0) then      
            write(iou_energies,'(1X, 1A8, 2X, 1A8, 2X, 1A37)') 'Time[ps],', & 
                & 'TimeStep[ms],', 'Kinetic energy[eV] of swapped atoms.'
            write(iou_energies,*)
        end if
        
        write(iou_energies,'(1X, 1F7.2, 1I7, 2F15.8)') t,  k, & 
            & convert(KEi(h), 'J','eV'), &
            & convert(KEi(c), 'J', 'eV')
        
    end subroutine writeStepEnergies
    
    subroutine writeBasicInfo ()
        real(wp):: start, finish
        integer::i
        
        call cpu_time(start)
        call buildSystem(lattice_const,latM,T0)
    	call doBox()
    	call cpu_time(finish)
        write(*,'(/, 1X,1A12, T40, 1A1, 2(1F4.1,", "), 1F5.1, 1A1)')'Box size[A]:', &
            &  '[',[(convert(box(i), 'm','A'),i=1,3)],']'
        write(*,'(1X, 1A16, T40, 1I5)') 'Number of steps:', N_steps
        write(*,'(1X, 1A17, T40, 1I5)') 'Skip_swap factor:', skip_swap
        write(*,'(1X, 1A16, T40, 1I5)') 'Number of atoms:', size(atoms)
        write(*,'(1X, 1A31, T40, 1I5)') 'Average number of atoms/region:', &
            & size(regionList(0.0_wp, real(latM(3)*lattice_const/N_slabs, wp)))
        write(*,'(1X, 1A25, T40, 1F4.1, 1A5)')'Time to build the system:', &
            & (finish-start), ' sec.'
        write(*,'(1X, 1A33, T40, 1A50)') 'Estimated time to run simulation:', & 
            writeUsedTime(estimateTime(N_steps))

        if (allocated(atoms)) deallocate(atoms)
        if (allocated(regions)) deallocate(regions)
        if (allocated(types)) deallocate(types)
            
    contains
        
        function estimateTime(N) result (o)
            integer, intent(in)::N
            real(wp):: o
            integer::k, iou_test, h, c
            
            open(file='mark1.test',newunit=iou_test) 
            
            call cpu_time(start)
            do k= 0, 10
                call rnem(k, h, c)
                if(mod(k,skip_dump)==0)     call writeStepXYZ(iou_test)
                if(mod(k,skip_neighbor)==0) call updateAllNeighbors()
                call velocityVerlet(dt)
                call doBox()
            end do
            call cpu_time(finish)
            
            o = (finish-start)*(N/11)
            
            close(iou_test)
            
        end function estimateTime
                       
    end subroutine writeBasicInfo
    
    function writeUsedTime(t) result (o)
        real(wp), intent(in)::t
        character(len=100)::o, secs, mins, hrs, days
        
        secs = '(1F5.2,1A5)'
        mins = '(2(1I2, 1A6))'
        hrs =  '(3(1I2, 1A6))'
        days = '(4(1I2, 1A6))'
        
        if (t <= 60) write(o,secs) t, ' sec.'
        if (t > 60 .and. t <= 3600) then 
            write(o,mins) nint(t/60), ' min, ', nint(mod(t,60.0_wp)), ' sec.'
        end if
        if (t > 3600 .and. t <= 86400) then 
            write(o,hrs) nint(t/60/60),   ' hrs, ', & 
                & nint(mod(t,60.0_wp)),   ' min, ', &
                & nint(mod(t/60,60.0_wp)),' sec.'
        end if
        
        if (t > 86400) then 
            write(o,days) nint(t/60/60/24),  ' days, ', & 
                & nint(mod(t,24.0_wp)),      ' hrs, ',  &
                & nint(mod(t/24,60.0_wp)),   ' min,',   & 
                & nint(mod(t/60/24,60.0_wp)),' sec.'       
        end if
        
    end function writeUsedTime

end module output_mod
