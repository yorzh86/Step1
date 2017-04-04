module system_mod
	use kinds_mod
	use utilities_mod
	use units_mod
	use settings_mod
	use autodiff_mod
	implicit none
	
	!=========!
	!= Types =!
	!=========!
	
	type::type_t
		!real(wp)::m
		type(ad_t)::m
			!! Atomic mass
		character(2)::atom_name
			!! Atom symbol (use standard periodic tables)
	end type
	
	type::atom_t
		type(ad_t),dimension(3)::r
			!! Atomic position
		type(ad_t),dimension(3)::v
			!! Atomic velocity
		type(ad_t),dimension(3)::a
			!! Atomic acceleration
		type(ad_t),dimension(3)::f
			!! Atomic force
		integer::atom_id
			!! Atom id
		integer::t
			!! Atom type
		integer,dimension(:),allocatable::neighbors
			!! Nearest neighbors
	end type
	
	type:: region_t
		type(ad_t),dimension(:),allocatable::temps
		type(ad_t),dimension(:),allocatable::Kenergies
		type(ad_t)::zl, zh
	end type
	
	
	!===============!
	!= Module Data =!
	!===============!
	
	logical::enableLennardJones = .false.
		!! Run Lennard-Jones potential
	
	type(type_t),dimension(:),allocatable::types
		!! Types for all atoms in system
	type(atom_t),dimension(:),allocatable::atoms
		!! All atoms in system
	type(region_t),dimension(:),allocatable::regions
	
	type(ad_t),dimension(:),allocatable::Totenergies
		!! array for total energy of the system at a timestep
	
	type(ad_t)::Teta
		!! Thermostat DOF
	type(ad_t)::Pepsilon
		!! Barostat DOF
	type(ad_t),dimension(3)::box
		!! Bounds of the simulation box
	
	integer,dimension(:),allocatable::listofRegions
	public::listofRegions
	
	integer::ts
		!! Time step counter
	type(ad_t)::t
		!! System time
	
	public:: KE
	public:: PE
	public:: E
	
contains

	subroutine buildSystem(a,N,Ti)
		implicit none
		type(ad_t), intent(in)::a
			!! Lattice constant
		integer,dimension(3),intent(in)::N
			!! Number of unit cells
		type(ad_t),intent(in)::Ti
			!! Initial temperature
		real(wp), dimension(3,4), parameter::rcell= &
		reshape([0.0_wp, 0.0_wp, 0.0_wp, &
				 0.5_wp, 0.5_wp, 0.0_wp, &
				 0.0_wp, 0.5_wp, 0.5_wp, &
				 0.5_wp, 0.0_wp, 0.5_wp  ], [3,4])
		integer::i,j,k,l,idx, ns
		
		Teta = 0.0_wp
		Pepsilon = 0.0_wp
		
		box = a*real(N,wp)
		ns = nint(real(N_steps/skip_swap,wp))
				
		allocate(types(1))
		allocate(atoms(size(rcell,2)*product(N)))
		allocate(regions(N_slabs))
		do i=1, N_slabs
			allocate(regions(i)%temps(0:N_steps))
		end do
		do i=1, 2
			allocate(regions(i)%Kenergies(0:N_steps/skip_swap))
		end do
		
		allocate(Totenergies(0:N_steps))
		
		!! Attention: units changed
		types%m = convert(39.948_wp,'u','u')
		
		types%atom_name = 'Ar'
		atoms(:)%t = 1
		
		idx = 1
		do k=1,N(3)
			do j=1,N(2)
				do i=1,N(1)
					do l=1,size(rcell,2)
						atoms(idx)%r = a*(real([i,j,k]-1,wp)+rcell(1:3,l))
						idx = idx+1
					end do
				end do
			end do
		end do
		
		do k=1,size(atoms)
			!! Create random direction of velocities
			call random_number(atoms(k)%v%x)  !must be %x
			atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			do while(norm2(real(atoms(k)%v))>1.0_wp .and. norm2(real(atoms(k)%v))<0.1_wp)
				call random_number(atoms(k)%v%x)
				atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			end do
			atoms(k)%v = atoms(k)%v/norm2(atoms(k)%v)
			!! Set velocity magnitude
			atoms(k)%v = atoms(k)%v*sqrt(2.0_wp*kB*Ti/types(atoms(k)%t)%m)*abs(randomNormal()+1.0_wp)
			!! velocity = self. * sqrt (eV/K * K / gram/mole) = self. * sqrt(eV/ g/mole)
			!! A/ps = eV/ u
		end do
		forall(k=1:3) atoms(:)%v(k) = (atoms(:)%v(k)-sum(atoms(:)%v(k))/real(size(atoms),wp))
		
		call updateAllLists()
		
		do k=1,size(atoms)
			atoms(k)%atom_id = k
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
			atoms(k)%f = -delV(k)
		end do
		
		ts = 0
		t  = 0.0_wp

	end subroutine buildSystem

	subroutine writeLammpsData(fn)
		character(*),intent(in)::fn
		
		integer::i,k,iou
		real(wp)::E0,S0

		
		E0 = lj%coeffs(1)
		S0 = lj%coeffs(2)
		
		open(file=fn,newunit=iou)
		
		write(iou,'(1A)') '# Input geometry for lammps'
		write(iou,'(1I9,1X,1A)') size(atoms),'atoms'
		write(iou,'(1I9,1X,1A)') size(types),'atom types'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,convert(box(1),'A','A'),'xlo','xhi'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,convert(box(2),'A','A'),'ylo','yhi'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,convert(box(3),'A','A'),'zlo','zhi'
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Masses'
		write(iou,'(1A)') ''
		do k=1,size(types)
			write(iou,'(1I4,1X,1F13.6)') k,convert(types(k)%m,'u','u')
		end do
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'PairIJ Coeffs # lj/cut'
		write(iou,'(1A)') ''
		write(iou,'(2I3,3F13.6)') 1, 1 , convert(E0,'uA2/ps2','eV'), convert(S0,'A','A'), convert(lj%cutoff,'A','A')
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Atoms'
		write(iou,'(1A)') ''
		do k=1,size(atoms)
			write(iou,'(1I9,1X,1I2,1X,3E25.15)') k,atoms(k)%t,[( convert(atoms(k)%r(i),'A','A') , i=1,3 )]
		end do
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Velocities'
		write(iou,'(1A)') ''
		do k=1,size(atoms)
			write(iou,'(1I9,1X,1X,3E25.15)') k,[( convert(atoms(k)%v(i),'A/ps','A/ps') , i=1,3 )]
		end do
		close(iou)
	end subroutine writeLammpsData

	function V(i) result(o)
	!! fn of potential
	!! o - should be energy [eV]
		integer,intent(in)::i
		type(ad_t)::o
		
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
		
	contains
	
		subroutine doLennardJones(o)
			type(ad_t),intent(out)::o
			type(ad_t),dimension(3)::d
			type(ad_t)::l,E0,S0
			integer::j,aj
			
			E0 = lj%coeffs(1)
			S0 = lj%coeffs(2)
			!! E0 - energy [eV]
			!! S0 - distance [A]
			
			do j=1,size(atoms(i)%neighbors)
				aj = atoms(i)%neighbors(j)
				!! aj - integer, number of neighbors
				d = deltaR(atoms(aj),atoms(i))
				!! d - distance [A]
				if( norm2(real(d))>real(lj%cutoff) ) cycle
				l = S0/norm2(d)
				!! l - unitless distance/ distance
				o = o+( 0.5_wp ) * ( 4.0_wp*E0*l**6*(l**6-1.0_wp) )
					!! Only include half of the bond energy for each atom in the bond
					!! o - energy [eV]
			end do
		end subroutine doLennardJones
	
	end function V

	function delV(i) result(o)
		!! Total force on atom [eV/A]
		integer,intent(in)::i
		type(ad_t),dimension(3)::o
		type(ad_t),dimension(3)::r
		integer::j,aj
		
		o = 0.0_wp
		do j=1,size(atoms(i)%neighbors)
			aj = atoms(i)%neighbors(j)
			!! aj - integer, number of neighbors
			r  = deltaR(atoms(i),atoms(aj))
			!! r - distance [A]
			if( norm2(real(r))>real(lj%cutoff) ) cycle
			o = o+delVij(i,aj,r)
			!! o = eV/A
			
		end do
	end function delV
	
	function delVij(i,j,d) result (o)
		!! Force between two atoms [eV/A]
		integer,intent(in)::i,j
		type(ad_t),dimension(3),intent(in)::d
		type(ad_t),dimension(3)::o
		
		o = 0.0_wp
		if(i==j) return
		
		if(enableLennardJones) call doLennardJones(o)
		
	contains
		
		subroutine doLennardJones(o)
		!! o - should be force between atoms for LJ pot [eV/A]
			type(ad_t),dimension(3),intent(out)::o
			type(ad_t)::l,E0,S0
			
			E0 = lj%coeffs(1)
			S0 = lj%coeffs(2)
			!! E0 - energy [eV]
			!! S0 - distance [A]
			
			l = S0/norm2(d)
			! l - unitless (distance / distance)
			o = o+24.0_wp*E0/sum(d*d)*(l**6)*(1.0_wp-2.0_wp*l**6)*d
			!! 24* eV / A^2 * (unitless^6) * (1 - 2 * unitless^6)* [A]
			!! eV/A
			
		end subroutine doLennardJones
		
	end function delVij

	function deltaR(a1,a2) result(o)
	!! distance between two atoms [A]
		type(atom_t),intent(in)::a1,a2
		type(ad_t),dimension(3)::o
		type(ad_t),dimension(3)::d
		d = a1%r-a2%r
		!! d - distance between two atoms (array of [A])
		o = d-box*real(nint(real(d/box)),wp)
		!! o - distance between atoms and box [A]
	end function deltaR

	function temperature() result(o)
	!! temperature [K]
		type(ad_t)::o
		type(ad_t),dimension(3)::vBulk
		type(ad_t)::SKE
		integer::k
		
		forall(k=1:3) vBulk(k) = sum(atoms%v(k))/real(size(atoms),wp)
		
		SKE = 0.0_wp
		!! sum of KE
		do k=1,size(atoms)
			SKE = SKE+KEi(k,vBulk)
		end do
		o = (2.0_wp*SKE)/(3.0_wp*real(size(atoms),wp)*kB)
		!! eV / (integer * eV/K) = [K]
	end function temperature

	function E() result(o)
	!! sum of energies of system [eV]
		type(ad_t)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+Ei(k)
		end do
	end function E

	function Ei(i) result(o)
	!! sum of energies of atom [eV]
		integer,intent(in)::i
		type(ad_t)::o
		
		o = KEi(i) + PEi(i) ! Do I need to do ABS?
	end function Ei

	function KE() result (o) 
	!! sum of kinetic energies of system [eV]
		type(ad_t)::o
		!real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+KEi(k)
		end do
	end function KE
	
	function KEi(i,vBulk) result (o)
	!! kinetic energy of atom [eV]
		integer,intent(in)::i
		type(ad_t),dimension(3),intent(in),optional::vBulk
		type(ad_t)::o
		type(ad_t),dimension(3)::v0
		
		v0 = 0.0_wp
		if(present(vBulk)) v0 = vBulk

		o = 0.5_wp*types(atoms(i)%t)%m*sum((atoms(i)%v-v0)**2)
		!! o = 0.5 * gram/mole * ([A]/[ps])^2 = eV ??!!!!!!! ATTENTION
		!! 1 eV = 1.602-19 [J]
		!! for 1 atom, 1m, 1m/sec, mass 39.948: 
		!! KE[SI] = 3.3168E-23
		!! KE[metal] = 1.9974E-23   
		!! 1.6605?
		
	end function KEi
	
	function PE() result (o)
	!! sum of potential energies of system [eV]
		type(ad_t)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+PEi(k)
		end do
	end function PE

	function PEi(i) result (o)
	!! potential energy of atom
		integer,intent(in)::i
		type(ad_t)::o
		
		o = V(i)
	end function PEi

	function Si(i) result(o)
		!! we dont even use it...
		integer,intent(in)::i
		type(ad_t),dimension(3,3)::o
		type(ad_t),dimension(3)::v,r,F
		integer::j,aj
		
		v = atoms(i)%v
		o = -types(atoms(i)%t)%m*dyadic(v,v)
		
		do j=1,size(atoms(i)%neighbors)
			aj = atoms(i)%neighbors(j)
			r  = deltaR(atoms(i),atoms(aj))
			if( norm2(real(r))>real(lj%cutoff) ) cycle
			F = delVij(i,aj,r)
			o = o-0.5_wp*(dyadic(r,F)+dyadic(F,r))
		end do
		
	end function Si

	subroutine updateAllLists()
		integer::k
		!real(wp)::abc
		!abc = real(latM(3)*lattice_const/N_slabs, wp)
		
		do k=1,size(atoms)
			call updateNeighbors(k)
		end do
!		do j=1, N_slabs
!			listofRegions = regionList(j*abc - abc, j*abc)
!		end do
	end subroutine updateAllLists

	subroutine updateNeighbors(i)
		integer,intent(in)::i
		type(ad_t),dimension(3)::r
		
		integer::k
		
		atoms(i)%neighbors = [integer::]
		do k=1,size(atoms)
			if(i==k) cycle
			r  = deltaR(atoms(i),atoms(k))
			if( norm2(real(r))>real(lj%cutoff+lj%skin) ) cycle
			atoms(i)%neighbors = [atoms(i)%neighbors,k]
		end do
	end subroutine updateNeighbors

	function averageNeighbors() result(o)
		type(ad_t)::o
		integer::k
		
		o = sum([( size(atoms(k)%neighbors) , k=1,size(atoms) )]) / real(size(atoms),wp)/2.0_wp
	end function averageNeighbors

	function virial() result(o)
		type(ad_t)::o
		type(ad_t),dimension(3)::F,r
		integer::i,j,aj
		
		o = 0.0_wp
		do i=1,size(atoms)
			do j=1,size(atoms(i)%neighbors)
				aj = atoms(i)%neighbors(j)
				r  = deltaR(atoms(i),atoms(aj))
				!! r - distance [A]
				if( norm2(real(r))>real(lj%cutoff) ) cycle
				F = delVij(i,aj,r)
				o = o+0.5_wp*sum(F*r)
				!! o  = self. + 0.5* (eV/A * A) = eV
			end do
		end do
	end function virial

	function pressure() result(o)
		type(ad_t)::o
		
		o = (real(size(atoms),wp)*kB*temperature()-virial()/3.0_wp)/(box(1)*box(2)*box(3))
		!! o = (atoms number * eV/K * K - (eV/A * A)/3)/A^3 = eV/A^3 = eV/A  / A^2 = force/area
	end function pressure

	function regionList(zl,zh) result(o)
		!! Stores an array of atoms(i)%atoms_id that belong to a region
		real(wp),intent(in)::zl,zh
		
		integer,dimension(:),allocatable::o
		
		integer::k
		
		o = pack( [( k, k=1,size(atoms))] , real(atoms%r(3))>=zl .and. real(atoms%r(3))<zh )
	end function regionList
	
	function listTemp(l) result(o)
		!! Calculates the average temperature of atoms in list L
		integer,dimension(:),intent(in)::l
		type(ad_t)::o
		
		
		integer::k,ai
		
		o = 0.0_wp
		do k=1,size(l)
			ai = l(k)
			o = o + types(atoms(ai)%t)%m*sum(atoms(ai)%v**2)
		end do
		
		o = o/(3.0_wp*kB*real(size(l),wp))
	end function listTemp
	
	function selectHot(l) result(o)
		!! Finds the hottest atom from the list of L
		integer,dimension(:),intent(in)::l
		integer::o
		
		real(wp),dimension(:),allocatable::T
		integer::k,ai
		
		allocate(T(size(l)))
		do k=1,size(l)
			ai = l(k)
			T(k) = 0.5_wp*types(atoms(ai)%t)%m*sum(atoms(ai)%v**2)/(3.0_wp*kB)
		end do
		
		o = l(maxloc(T,1))
	end function selectHot
	
	function selectCold(l) result(o)
		!! Finds the coldest atomf from list L
		integer,dimension(:),intent(in)::l
		integer::o
		
		real(wp),dimension(:),allocatable::T
		integer::k,ai
		
		allocate(T(size(l)))
		do k=1,size(l)
			ai = l(k)
			T(k) = 0.5_wp*types(atoms(ai)%t)%m*sum(atoms(ai)%v**2)/(3.0_wp*kB)
		end do
		
		o = l(minloc(T,1))
	end function selectCold
	
end module system_mod
