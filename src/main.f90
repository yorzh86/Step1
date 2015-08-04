program main_prg
	use kinds_mod
	use system_mod
	use integrate_mod
	use output_mod
	implicit none
	
	real(wp)::T0 = 40.0_wp
	real(wp)::dt
	integer::iou,k
	
	open(file='out.xyz',newunit=iou)
	
	enableLennardJones = .true.
	call setThermostat(.true.,T0,100.0_wp)
	dt = 0.1_wp
	
	call buildSystem(5.260_wp,15,T0)
	!(lattice parameter, box edge, temperature)
	call doBox()
	call writeStepXYZ(iou)
	
	do k=1,150
		call velocityVerlet(dt)
		call doBox()
		if(k==10000) call setThermostat(.false.)
		
		if(mod(k,1)==0) then
			call writeStepXYZ(iou)
			write(*,*) k,temperature(),PE(),KE(),heatflux()
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
