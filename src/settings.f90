module settings_mod
	use kinds_mod
	use units_mod
	use autodiff_mod
	implicit none
	
	!=======================!
	!= Universal Constants =!
	!=======================!
	
	real(wp),parameter::kB = 0.83144618021123218  !g/mole*A^2/ps^2/K
	!!kg to g/mole = 6.022140857E+26
	!!m2 to A2= 1E+20
	!!s2 to ps2 = 1E+24
	!!J = N * m; N = kg * m/s2; J = kg * m2/s2
	!!8.617342371E-5_wp ![eV/K]
	!!1.3806488E-23_wp  ![J/K] 
	!! Boltzmann constant
	
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
		
		E0 = kB*convert(125.7_wp,'K','K')
		S0 = convert(3.345_wp, 'A','A')
		lj%cutoff = 2.5_wp*S0
		lj%skin = 0.5_wp*S0
		
		lj%coeffs = [diff(E0,1),diff(S0,2)]
		
		!= Simulation =!
		N_steps       = 100000
		N_slabs       = 10
		skip_swap     = 500
		skip_thermo   = 100
		skip_dump     = 5
		skip_neighbor = 200
		
		lattice_const = convert(5.40_wp, 'A', 'A')
		latM = [5,5,5]
		
		T0 = convert(45.0_wp,'K','K')
		P0 = convert(1.0_wp,'bar','bar')
		dt = convert(0.01_wp,'ps','ps')!10fs = 0.01 ps
		
		!= Thermostat =!
		thermostat%tau = convert(100.0_wp*dt,'s','ps')
		
		!= Barostat =!
		barostat%tau = convert(1000.0_wp*dt,'s','ps')
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

		write(iou,'(1A,1EN25.5)') 'variable dt            equal ',convert(dt,'ps','ps')
		write(iou,'(1A,1EN25.5)') 'variable tau_T         equal ',convert(thermostat%tau,'ps','ps')
		write(iou,'(1A,1EN25.5)') 'variable tau_P         equal ',convert(barostat%tau,'ps','ps')
		write(iou,'(1A,1EN25.5)') 'variable T0            equal ',convert(T0,'K','K')
		write(iou,'(1A,1EN25.5)') 'variable P0            equal ',convert(P0,'bar','bar')
		
		write(iou,'(1A,1EN25.5)') 'variable cutoff        equal ',convert(lj%cutoff,'A','A')
		write(iou,'(1A,1EN25.5)') 'variable skin          equal ',convert(lj%skin,'A','A')
		
		write(iou,'(1A,1EN25.5)') 'variable box_length    equal ',real(latM(3),wp)*convert(lattice_const,'A','A')
		
		close(iou)
	end subroutine writeLammpsVars

end module settings_mod
