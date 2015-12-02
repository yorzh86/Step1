module settings_mod
	use kinds_mod
	use units_mod
	implicit none
	
	!=======================!
	!= Universal Constants =!
	!=======================!
	
	real(wp),parameter::kB = 1.3806488E-23_wp
		!! Boltzmann constant in SI units
	
	!=========!
	!= Types =!
	!=========!
	
	type::potential_t
		real(wp)::cutoff
			!! Cutoff radius
		real(wp)::skin
			!! Skin beyond cutoff for neighbor lists
		real(wp),dimension(:),allocatable::coeffs
			!! Potential coefficients
	end type
	
	type::nosehoover_t
		real(wp)::eta = 0.0_wp
			!!
		real(wp)::set = 10.0_wp
			!!
		real(wp)::tau = 100.0_wp
			!!
	end type
	
	!=============!
	!= Variables =!
	!=============!
	
	type(potential_t)::lj ! coeffs=[E0,S0]
		!! Lennard-Jones potential data
	type(nosehoover_t)::thermostat
	
	integer::N_steps
		!! Number of total simulation steps
	integer::skip_thermo
		!! Number of steps between thermo reports
	integer::skip_dump
		!! Number of steps between xyz file dumps
	integer::skip_neighbor
		!! Number of steps between neighbor list rebuilds
	
	real(wp)::T0
		!! Initial temperature [K]
	real(wp)::dt
		!! Timestep [seconds]
	
contains

	subroutine initialize()
		real(wp)::E0,S0
		
		!= Potential =!
		
		E0 = kB*convert(125.7_wp,'K','K')
		S0 = convert(3.345_wp,'A','m')
		lj%cutoff = 2.0_wp*S0
		lj%skin = convert(1.0_wp,'A','m')
		
		lj%coeffs = [E0,S0]
		
		!= Simulation =!
		N_steps       = 3000
		skip_thermo   = 100
		skip_dump     = 100
		skip_neighbor = 20
		
		T0 = convert(40.0_wp,'K','K')
		dt = convert(10.0_wp,'fs','s')
		
		!= Thermostat =!
		thermostat%eta = 0.0_wp
		thermostat%tau = convert(100.0_wp*dt,'s','ps')
	end subroutine initialize

	subroutine writeLammpsVars(fn)
		character(*),intent(in)::fn
		integer::iou
		
		open(file=fn,newunit=iou)
		
		write(iou,'(1A,1I10)')   'variable N_steps       equal ',N_steps
		write(iou,'(1A,1I10)')   'variable skip_thermo   equal ',skip_thermo
		write(iou,'(1A,1I10)')   'variable skip_dump     equal ',skip_dump
		write(iou,'(1A,1I10)')   'variable skip_neighbor equal ',skip_neighbor
		
		write(iou,'(1A,1E25.5)') 'variable dt            equal ',convert(dt,'s','ps')
		write(iou,'(1A,1E25.5)') 'variable dt_nh         equal ',convert(thermostat%tau,'s','ps')
		write(iou,'(1A,1E25.5)') 'variable T0            equal ',convert(T0,'K','K')
		
		write(iou,'(1A,1E25.5)') 'variable cutoff        equal ',convert(lj%cutoff,'m','A')
		write(iou,'(1A,1E25.5)') 'variable skin          equal ',convert(lj%skin,'m','A')
		
		write(iou,'(1A)')        'variable t             equal step*${dt}'
		write(iou,'(1A)')        'variable T             equal temp'
		
		close(iou)
	end subroutine writeLammpsVars

end module settings_mod