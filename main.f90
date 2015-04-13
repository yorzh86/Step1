program main_prg
	use kinds_mod
	use system_mod
	use integrate_mod
	use output_mod
	implicit none
	
	real(wp)::T0 = 20.0_wppymo	
	real(wp)::dt
	integer::iou,k
	
	open(file='out.xyz',newunit=iou)
	
	enableLennardJones = .true.
	call setThermostat(.true.,T0,100.0_wp)
	dt = 1.0_wp
	
	!subroutine buildSystem(lattice parameter,box edge,temperature)
	call buildSystem(3.1_wp,10,T0)
	
	do k=1,50000
		call velocityVerlet(dt)
		call doBox()
		
		if(k==10000) then
			call setThermostat(.false.)
		end if
		
		if(mod(k,100)==0) then
			call writeStepXYZ(iou)
			write(*,*) k,temperature()
		end if
	end do
	
	close(iou)
	
contains



end program main_prg 

!================ metal units ==================
!mass              = grams/mole
!distance          = Angstroms
!time              = picoseconds
!energy            = eV
!velocity          = Angstroms/picosecond
!force             = eV/Angstrom
!torque            = eV
!temperature       = Kelvin
!pressure          = bars
!charge            = multiple of electron charge
!dipole            = charge*Angstroms
!density           = gram/cm^dim
!dynamic viscosity = Poise
!electric field    = volts/Angstrom
