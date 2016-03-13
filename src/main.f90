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
    integer::iou_lammps=1
        !! I/O unit to read lammps dump file
    integer::iou_energies
        !! I/O unit to write energies before swap
    integer::iou_temps
        !! I/O unit to write temperatures
        
    call setupSim()
    call runSim()
    call endSim()
    
contains

    subroutine setupSim
        open(file='mark1.xyz',newunit=iou_xyz)        
        open(file='mark1.temps',newunit=iou_temps)
        open(file='mark1.energies',newunit=iou_energies)
                        
        call initialize_parameters()
                    
        enableLennardJones = .true.
        call setThermostat(.false.,T0,10.0_wp*dt)
        call setBarostat(.false.,P0, 5.0E10_wp*dt)
        call buildSystem(lattice_const,latM,T0)
        
        call doBox()
        call writeLammpsData('Ar.data')
        call writeLammpsVars('Ar.vars')
    end subroutine setupSim
        
    subroutine runSim
        integer::k
        integer,dimension(:), allocatable::l
        
        write(*,'(1X, 1A12, 3F10.2)') 'Box size[A]:', box*1E10
        write(*,'(1X, 1A16, 1I6)') 'Number of steps:', N_steps
		write(*,'(1X, 1A16, 1I6)') 'Number of atoms:', size(atoms)
		write(*,'(1X, 1A31, 1I6)') 'Average number of atoms/region:', &
			& size(regionList(0.0_wp, real(latM(3)*lattice_const/N_slabs, wp)))
		
		write(*,'(1X, 1A37, 1F7.1, 1A4)') 'Estimated time to execute simulation:', & 
			& N_steps/10/60.0, 'min' 
                
        do k=0, N_steps
            call writeStepThermo(k, iou_temps, iou_energies)
            if(mod(k,skip_dump)==0)     call writeStepXYZ(iou_xyz)
            if(mod(k,skip_neighbor)==0) call updateAllNeighbors()              
            call velocityVerlet(dt)
            call doBox()
        end do
        write(*,*)
        write(*,*)"FINISHED!!!"
        
    end subroutine runSim

    subroutine endSim
        close(iou_xyz)
        close(iou_thermo)
        close(iou_temps)
        close(iou_energies)
    end subroutine endSim
    
end program main_prg 
