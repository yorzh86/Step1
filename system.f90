module system_mod
	use kinds_mod
	implicit none
	
	!==============!
	!= Parameters =!
	!==============!
	

	real(wp),parameter::kB = 8.617332478E-5_wp
	real(wp),parameter::E0 = 1.0298490416E-2_wp
	real(wp),parameter::S0 = 3.4_wp 
	
	!=========!
	!= Types =!
	!=========!
	
	type::type_t
		real(wp)::m
		character(2)::atom_name
	end type
	
	type::atom_t
		real(wp),dimension(2)::r
		real(wp),dimension(2)::v
		real(wp),dimension(2)::a
		integer::t
	end type
	
	!===============!
	!= Module Data =!
	!===============!
	
	logical::enableLennardJones = .false.
	
	type(type_t),dimension(:),allocatable::types
	type(atom_t),dimension(:),allocatable::atoms
	
	real(wp),dimension(2)::box
	
	integer::ts
	real(wp)::t
	
contains

	subroutine buildSystem(a,N,Ti)
	! creates a system of atoms using lattice par, box edge length and
	! temperature
		real(wp),intent(in)::a
		integer,intent(in)::N
		real(wp),intent(in)::Ti
		
		integer::i,j,k
		
		box = a*real([N,N],wp)*[1.0_wp,1.0_wp/sqrt(2.0_wp)]
		
		allocate(types(1))
		allocate(atoms(N**2))
		
		types%m = 39.9_wp
		types%atom_name = 'Ar'
		atoms(:)%t = 1
		
		forall(i=1:N,j=1:N)
		! sets initial position and acceleration of atoms
			atoms(i+N*(j-1))%r = a*(real([i,j],wp)*[1.0_wp,1.0_wp/sqrt(2.0_wp)]+[mod(j,2)/2.0_wp,0.0_wp]) &
			& -a/2.0_wp+a/2.0_wp*N
			atoms(i+N*(j-1))%a = -delV(i+N*(j-1))/types(atoms(i+N*(j-1))%t)%m
		end forall
		
		do k=1,N**2
		! applies normal distribution for initial velocities
			call random_number(atoms(k)%v)
			atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			do while(norm2(atoms(k)%v)>1.0_wp .and. norm2(atoms(k)%v)<0.1_wp)
				call random_number(atoms(k)%v)
				atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			end do
			atoms(k)%v = atoms(k)%v/norm2(atoms(k)%v)
			atoms(k)%v = atoms(k)%v*sqrt(2.0_wp*kB*Ti/types(atoms(k)%t)%m)*abs(randomNormal()+1.0_wp)
		end do
		
		forall(k=1:2) atoms(:)%v(k) = atoms(:)%v(k)-sum(atoms(:)%v(k))/real(N*N,wp)
		
		ts = 0
		t  = 0.0_wp
	end subroutine buildSystem

	pure function V(i) result(o)
	! helps to calculate Epot at different potentials
		integer,intent(in)::i
		real(wp)::o
		
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
		
	contains
	
		pure subroutine doLennardJones(o)
		
			real(wp),intent(inout)::o
			
			real(wp),dimension(2)::d
			real(wp)::l
			integer::k
			
			do k=1,size(atoms)
				if(k==i) cycle
				d = deltaR(atoms(k),atoms(i))
				l = S0/norm2(d)
				o = 4*l**12-l**6
			end do
		end subroutine doLennardJones
	
	end function V

	pure function delV(i) result(o)
	! calculates forces between atoms at employed potentials
		integer,intent(in)::i
		real(wp),dimension(2)::o
		integer::j
		
		o = 0.0_wp
		do j=1,size(atoms)
			o = o+delVij(i,j)
		end do
		
	end function delV

	function delVij(i,j)
		integer,intent(in)::i,j
		real(wp),dimension(2)::o
		
		o = 0.0_wp
		if(i==j) return
		
		if(enableLennardJones) call doLennardJones(o)
		
	contains
		
		pure subroutine doLennardJones(o)
			real(wp),dimension(2),intent(inout)::o
			
			real(wp),dimension(2)::d
			real(wp)::l
			
			d = deltaR(atoms(i),atoms(j))
			l = S0/norm2(d)
			o = o+24.0_wp*E0/sum(d*d)*(2.0_wp*l**12-l**6)*d
		end subroutine doLennardJones
		
	end function delVij

	pure function deltaR(a1,a2) result(o)
	! calculates the distance between atoms, applying 
	! minimum image convention (nint(d/box))
		type(atom_t),intent(in)::a1,a2
		real(wp),dimension(2)::o
		
		real(wp),dimension(2)::d
		
		d = a1%r-a2%r
		o = d-box*nint(d/box)
	end function deltaR

	pure function temperature() result(o)
	! calculates temperature of atoms
		real(wp)::o
		integer::k
		
		! Change to subtract off average velocity
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+0.5_wp*types(atoms(k)%t)%m*norm2(atoms(k)%v)**2
		end do
		o = o/(2.0_wp*real(size(atoms),wp)*kB)
	end function temperature

	pure function KE() result (o)
	! calculates Ekin 
		real(wp)::o
		integer::k
		o = 0.0_wp
		do k=1,size(atoms)
			o = o + 0.5_wp*types(atoms(k)%t)%m*norm2(atoms(k)%v)**2
		end do
	end function KE
	
	pure function KEi(i) result (o)
		integer,intent(in)::i
		real(wp)::o
		integer::k
		o = 0.5_wp*types(atoms(i)%t)%m*norm2(atoms(i)%v)**2
	end function KEi
	
	pure function PE() result (o)
	! calculates Epot
		real(wp)::o
		real(wp),dimension(2)::d
		real(wp)::l
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+V(k)
		end do
	end function PE

	pure function PEi(i) result (o)
	! why do we need to separate PEi and KEi functions?
	! why don't we calculate site energy in one function?
		integer,intent(in)::i
		real(wp),dimension(2)::o
		!real(wp)::o
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
	
	contains
		pure subroutine doLennardJones(o)
			real(wp),dimension(2),intent(inout)::o
			real(wp),dimension(2)::d
			real(wp)::l
			d = deltaR(atoms(i-1),atoms(i))
			l = norm2(d)
			o = 4*l**12-l**6
		end subroutine doLennardJones
	end function PEi
	
	subroutine heatfluxJ(i)
		integer,intent(in)::i
		real(wp),dimension(2)::o
		
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
		write(*,*) o
	contains
		pure subroutine doLennardJones(o)
			real(wp),dimension(2),intent(inout)::o
			real(wp),dimension(2)::d, fij, first_term, second_term
			real(wp)::l
			integer:: k,j
			
			d = deltaR(atoms(i+1),atoms(i))
			l = S0/norm2(d)
			fij = 24.0_wp*E0/sum(d*d)*(2.0_wp*l**12-l**6)*d
			
			first_term = real([0,0],wp)
			second_term = real([0,0],wp)
			
			do k=1, size(atoms)
				first_term = first_term +(KEi(i)+PEi(i))*atoms(i)%v
			end do
			
			do k=1, size(atoms)
				second_term = second_term + (fij*atoms(i+1)%v)*d
			end do
			
			o = (second_term + first_term)/box	

		end subroutine doLennardJones
	end subroutine heatfluxJ

end module system_mod
