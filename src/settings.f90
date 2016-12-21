module settings_mod
	use kinds_mod
	use units_mod
	use autodiff_mod
	implicit none
	
	!=======================!
	!= Universal Constants =!
	!=======================!
	
	real(wp),parameter::kB = 8.617342371E+5_wp!1.3806488E-23_wp [J/K]
		!! Boltzmann constant in metal units
	
	!=========!
	!= Types =!
	!=========!
	
	type::potential_t
		!real(wp)::cutoff
		type(ad_t)::cutoff
			!! Cutoff radius
		!real(wp)::skin
		type(ad_t)::skin
			!! Skin beyond cutoff for neighbor lists
		!type(ad_t),dimension(:),allocatable::coeffs
		type(ad_t),dimension(2)::coeffs
			!! Potential coefficients       
	end type
	
	type::nosehoover_t
		!real(wp)::set = 10.0_wp
		type(ad_t)::set
			!! Set value
		!real(wp)::tau = 100.0_wp
		type(ad_t)::tau
			!! Characteristic time
	end type
	
	!=============!
	!= Variables =!
	!=============!
	
	type(potential_t)::lj ! coeffs=[E0,S0]
		!! Lennard-Jones potential data
	type(nosehoover_t)::thermostat
		!! Thermostat settings
	type(nosehoover_t)::barostat
		!! Barostat settings
	
	integer::N_steps
		!! Number of total simulation steps
	integer:: N_slabs
		!! Number of domains perpendicular to z axis
	integer::skip_thermo
		!! Number of steps between thermo reports
	integer::skip_dump
		!! Number of steps between xyz file dumps
	integer::skip_neighbor
		!! Number of steps between neighbor list rebuilds
	integer::skip_swap
		!! Number of steps between swapping atoms velocities
		
	!real(wp)::lattice_const
	type(ad_t)::lattice_const
		!! Lattice constant [A]
	integer, dimension(3)::latM
		!! Lattice multiplier at buildSystem phase
	
	!real(wp)::T0
	type(ad_t)::T0
		!! Target temperature [K]
	!real(wp)::P0
	type(ad_t)::P0
		!! Target pressure [bar]
	!real(wp)::dt
	type(ad_t)::dt
		!! Timestep [picoseconds]
	
contains

	subroutine initialize_parameters()
		real(wp)::E0, S0
		
		!= Potential =!
		
		E0 = kB*125.7_wp
		S0 = 3.345_wp ![A]
		lj%cutoff = 2.5_wp*S0
		lj%skin = 0.5_wp*S0
		
		lj%coeffs = [diff(E0,1),diff(S0,2)]
		
		!= Simulation =!
		N_steps       = 30
		N_slabs       = 10
		skip_swap     = 5
		skip_thermo   = 100
		skip_dump     = 1
		skip_neighbor = 200
		
		lattice_const = 5.40_wp !not converted
		latM = [5,5,10]
		
		T0 = 45.0_wp !K
		P0 = 1.0_wp !bar
		dt = convert(10.0_wp,'fs','ps')
		
		!= Thermostat =!
		thermostat%tau = 100.0_wp*dt
		!thermostat%tau = convert(100.0_wp*dt,'s','ps')
		
		!= Barostat =!
		barostat%tau = 1000.0_wp*dt
		!barostat%tau = convert(1000.0_wp*dt,'s','ps')
	end subroutine initialize_parameters

	subroutine writeLammpsVars(fn)
		character(*),intent(in)::fn
		integer::iou
		
		open(file=fn,newunit=iou)
		
		write(iou,'(1A,1I10)')   'variable N_steps       equal ',N_steps
		write(iou,'(1A,1I10)')   'variable skip_thermo   equal ',skip_thermo
		write(iou,'(1A,1I10)')   'variable skip_dump     equal ',skip_dump
		write(iou,'(1A,1I10)')   'variable skip_neighbor equal ',skip_neighbor
		write(iou,'(1A,1I10)')   'variable N_slabs       equal ',N_slabs
		write(iou,'(1A,1I10)')   'variable skip_swap     equal ',skip_swap
		! all converts are cleared
		write(iou,'(1A,1EN25.5)') 'variable dt            equal ',dt
		write(iou,'(1A,1EN25.5)') 'variable tau_T         equal ',thermostat%tau
		write(iou,'(1A,1EN25.5)') 'variable tau_P         equal ',barostat%tau
		write(iou,'(1A,1EN25.5)') 'variable T0            equal ',T0
		write(iou,'(1A,1EN25.5)') 'variable P0            equal ',P0
		
		write(iou,'(1A,1EN25.5)') 'variable cutoff        equal ',lj%cutoff
		write(iou,'(1A,1EN25.5)') 'variable skin          equal ',lj%skin
		
		write(iou,'(1A,1EN25.5)') 'variable box_length    equal ',real(latM(3),wp)*lattice_const
		
		close(iou)
	end subroutine writeLammpsVars

end module settings_mod
